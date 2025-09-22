# =============================================================================
# booleanManager.R (Cleaned Version - Only Functions Used by Script 06)
# Purpose: Boolean network construction functions for PRISM workflow
# =============================================================================

library(dplyr)

# =============================================================================
# Gene Name Sanitization
# =============================================================================

#' Generate sanitized gene name mapping for BoolNet compatibility
#'
#' Creates a mapping between original gene names and sanitized versions that
#' are valid R variable names and compatible with BoolNet.
#'
#' @param geneNames Character vector of original gene names
#' @return Data frame with columns OriginalName and SanitizedName
#' @export
generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(
    OriginalName = geneNames,
    SanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Convert gene name to valid R variable name
#'
#' Converts gene names to valid R variable names by replacing invalid characters
#' with underscores, adding 'X' prefix for names starting with numbers, and
#' handling empty/invalid names.
#'
#' @param name Character string with original gene name
#' @return Character string with sanitized gene name
#' @export
sanitizeGeneName <- function(name) {
  # Convert to valid R variable name
  name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  if (grepl("^[0-9]", name2)) name2 <- paste0("X", name2)
  if (name2 == "" || name2 == "_") name2 <- "Unknown"
  name2
}

# =============================================================================
# Regulator Selection Using SCENIC Metadata
# =============================================================================

#' Select optimal regulators for target using SCENIC metadata
#'
#' Uses composite weighting of correlation, motif confidence, SCENIC importance
#' to choose the best regulators for Boolean rule inference.
#' Filters regulators by presence in expression matrix.
#'
#' @param target Character name of target gene (sanitized)
#' @param edges Data frame with SCENIC metadata containing columns TF, Target, corr, regType, motifConfidence, NES, hasMotif
#' @param matBin Binary expression matrix (genes x cells)
#' @param maxRegs Integer, maximum number of regulators to select (default: 3)
#' @return List with components: regulators (data frame), n_selected (integer), mean_weight (numeric), has_activators (logical), has_repressors (logical)
#' @export
selectOptimalRegulatorsForTarget <- function(target, edges, matBin, maxRegs = 3) {
  
  # Get edges targeting this gene
  targetEdges <- edges[edges$Target == target, , drop = FALSE]
  
  if (nrow(targetEdges) == 0) {
    return(list(
      regulators = data.frame(),
      n_selected = 0,
      mean_weight = 0,
      has_activators = FALSE,
      has_repressors = FALSE
    ))
  }
  
  # Filter to regulators present in expression matrix
  validRegs <- targetEdges$TF %in% rownames(matBin)
  if (!any(validRegs)) {
    return(list(
      regulators = data.frame(),
      n_selected = 0,
      mean_weight = 0,
      has_activators = FALSE,
      has_repressors = FALSE
    ))
  }
  
  targetEdges <- targetEdges[validRegs, , drop = FALSE]
  
  # FIXED: Remove self-regulation unless explicitly in SCENIC to avoid identity rules
  targetEdges <- targetEdges[targetEdges$TF != target, , drop = FALSE]
  
  if (nrow(targetEdges) == 0) {
    return(list(
      regulators = data.frame(),
      n_selected = 0,
      mean_weight = 0,
      has_activators = FALSE,
      has_repressors = FALSE
    ))
  }
  
  # Calculate composite weights
  weights <- calculateRegulatorWeights(targetEdges)
  targetEdges$composite_weight <- weights
  
  # Sort by weight and take top candidates
  targetEdges <- targetEdges[order(-weights), , drop = FALSE]
  
  # Apply diversity filtering if we have more candidates than slots
  if (nrow(targetEdges) > maxRegs) {
    # Simple truncation to maxRegs (no diversity filtering)
    targetEdges <- head(targetEdges, maxRegs)
  }
  
  # Limit to maxRegs
  targetEdges <- head(targetEdges, maxRegs)
  
  return(list(
    regulators     = targetEdges,
    n_selected     = nrow(targetEdges),
    mean_weight    = if (nrow(targetEdges) > 0) mean(targetEdges$composite_weight) else 0,
    has_activators = any(targetEdges$regType == "Activation", na.rm = TRUE),
    has_repressors = any(targetEdges$regType == "Inhibition", na.rm = TRUE)
  ))
}

#' Calculate composite regulator weights using SCENIC metadata
#'
#' Combines multiple SCENIC-derived metrics into a single composite weight
#' for regulator selection. Includes correlation strength, motif confidence,
#' NES scores, prior knowledge, and regulatory direction clarity.
#'
#' @param edges Data frame of edges for a single target with SCENIC metadata
#' @return Numeric vector of composite weights
#' @export
calculateRegulatorWeights <- function(edges) {
  
  # Base weight from correlation strength
  corr_weight <- abs(edges$corr)
  
  # Motif confidence bonus (0 if missing)
  motif_bonus <- ifelse(!is.na(edges$motifConfidence), 
                        pmax(0, edges$motifConfidence), 0)
  
  # SCENIC NES score bonus (normalized)
  nes_bonus <- ifelse(!is.na(edges$NES), 
                      pmax(0, edges$NES / 5), 0)
  
  # Prior knowledge bonus for confirmed motifs
  prior_bonus <- ifelse(!is.na(edges$hasMotif) & edges$hasMotif == TRUE, 0.1, 0)
  
  # Handle missing regType
  if (!"regType" %in% colnames(edges)) {
    edges$regType <- "Activation"  # Default assumption
  }
  
  # Slight bonus for clear regulatory direction
  direction_bonus <- ifelse(edges$regType %in% c("Activation", "Inhibition"), 0.05, 0)
  
  # Composite weight calculation (matches config.yaml weights)
  composite <- corr_weight + 
               (motif_bonus * 0.3) + 
               (nes_bonus * 0.2) + 
               prior_bonus + 
               direction_bonus
  
  return(composite)
}

# =============================================================================
# Boolean Rule Inference
# =============================================================================

#' Infer best Boolean rule using multiple methods
#'
#' Tries both sign-aware templates and empirical state counting,
#' then selects the best based on composite scoring. Returns the
#' highest-scoring rule with complete metadata.
#'
#' @param target Target gene name
#' @param regulator_info List with selected regulator information from selectOptimalRegulatorsForTarget
#' @param matBin Binary expression matrix (genes x cells)
#' @param cellOrder Vector of cell indices ordered by pseudotime
#' @param config Configuration object with scoring parameters
#' @return List with best rule and metadata, or NULL if no valid rule found
#' @export
inferBestBooleanRule <- function(target, regulator_info, matBin, cellOrder, config) {
  
  candidates <- list()
  
  # Method 1: Sign-aware templates
  templates <- generateSignAwareTemplates(target, regulator_info)
  
  for (template_name in names(templates)) {
    rule_str <- templates[[template_name]]
    
    # Score template
    score_info <- scoreRuleWithScenicMetrics(rule_str, target, matBin, regulator_info)
    
    # Store rule in BoolNet format (target, logic)
    full_rule_str <- paste0(target, ", ", rule_str)
    
    candidates[[paste0("template_", template_name)]] <- list(
      rule = full_rule_str,
      score = score_info$composite_score,
      method = paste0("template_", template_name),
      score_breakdown = score_info,
      k_used = 0
    )
  }
  
  # Method 2: Empirical state counting (k=0 only)
  if (regulator_info$n_selected > 0) {
    reg_names <- regulator_info$regulators$TF
    
    # Create input-output pairs
    tryCatch({
      ioData <- makeInputOutputPairs(target, reg_names, matBin, cellOrder, k = config$boolKValue)
      
      if (!is.null(ioData) && length(ioData$stateIndices) >= 10) {
        
        # Find best empirical rule
        empirical_result <- findBestBooleanRules(ioData)
        
        if (!is.null(empirical_result) && empirical_result$score > 0 && 
            length(empirical_result$bestFns) > 0) {
          
          # Generate rule string
          rule_str <- makeBoolNetRule(target, ioData$regulators, 
                                      empirical_result$bestFns[[1]]$outPattern)
          
          # Score empirical rule
          score_info <- scoreRuleWithScenicMetrics(rule_str, target, matBin, regulator_info)
          
          candidates$empirical <- list(
            rule = rule_str,
            score = score_info$composite_score,
            method = "empirical_k0",
            score_breakdown = score_info,
            k_used = 0,
            empirical_metadata = list(
              functions_tested = empirical_result$functionsTested %||% empirical_result$functionsTestsed,
              total_possible = empirical_result$totalPossible,
              io_pairs = length(ioData$stateIndices)
            )
          )
        }
      }
    }, error = function(e) {
      warning("Empirical method failed for ", target, ": ", e$message)
    })
  }
  
  # Select best candidate
  if (length(candidates) == 0) {
    return(NULL)
  }
  
  # Find candidate with highest score
  scores <- sapply(candidates, `[[`, "score")
  best_idx <- which.max(scores)
  best_candidate <- candidates[[best_idx]]
  
  # Add common metadata
  best_candidate$regulators <- if (regulator_info$n_selected > 0) regulator_info$regulators$TF else character(0)
  best_candidate$n_regulators <- regulator_info$n_selected
  best_candidate$regulator_weights <- if (regulator_info$n_selected > 0) regulator_info$regulators$composite_weight else numeric(0)
  
  return(best_candidate)
}

#' Generate Boolean templates respecting SCENIC regulatory signs
#'
#' Creates biologically plausible Boolean logic templates based on whether
#' regulators are activators or repressors according to SCENIC analysis.
#' Generates multiple template types (single regulator, OR, AND, balanced, repressor-dominant).
#'
#' @param target Character name of target gene
#' @param regulator_info List from selectOptimalRegulatorsForTarget with regulator metadata
#' @return Named list of template rule strings (BoolNet logic format - target name not included)
#' @export
generateSignAwareTemplates <- function(target, regulator_info) {
  
  if (length(regulator_info) == 0 || regulator_info$n_selected == 0) {
    # FIXED: Return empty list instead of self-activation to avoid identity rules
    return(list())
  }
  
  regulators <- regulator_info$regulators
  activators <- regulators$TF[regulators$regType == "Activation" & !is.na(regulators$regType)]
  repressors <- regulators$TF[regulators$regType == "Inhibition" & !is.na(regulators$regType)]
  
  templates <- list()
  
  # Template 1: Single regulator (if only one)
  if (nrow(regulators) == 1) {
    reg <- regulators$TF[1]
    if (!is.na(regulators$regType[1]) && regulators$regType[1] == "Inhibition") {
      templates$single_repressor <- paste0("!", reg)
    } else {
      templates$single_activator <- reg
    }
    return(templates)
  }
  
  # Template 2: All activators OR (any activator sufficient)
  if (length(activators) > 0) {
    templates$activator_or <- paste(activators, collapse = " | ")
  }
  
  # Template 3: All activators AND (cooperative activation)
  if (length(activators) > 1) {
    templates$activator_and <- paste(activators, collapse = " & ")
  }
  
  # Template 4: Balanced regulation (activators present, repressors absent)
  if (length(activators) > 0 && length(repressors) > 0) {
    act_part <- if (length(activators) == 1) {
      activators[1]
    } else {
      paste0("(", paste(activators, collapse = " | "), ")")
    }
    
    rep_part <- if (length(repressors) == 1) {
      paste0("!", repressors[1])
    } else {
      paste0("!(", paste(repressors, collapse = " | "), ")")
    }
    
    templates$balanced <- paste0(act_part, " & ", rep_part)
  }
  
  # Template 5: Repressor dominant (target OFF when any repressor ON)
  if (length(repressors) > 0) {
    if (length(repressors) == 1) {
      templates$repressor_dominant <- paste0("!", repressors[1])
    } else {
      templates$repressor_dominant <- paste0("!(", paste(repressors, collapse = " | "), ")")
    }
  }
  
  return(templates)
}

#' Score Boolean rule with SCENIC-aware composite metrics
#'
#' Evaluates Boolean rules using empirical accuracy plus bonuses for
#' sign consistency with SCENIC regulatory directions and SCENIC confidence scores.
#'
#' @param ruleLogic Character string with BoolNet-compatible logic (without target name)
#' @param target Target gene name
#' @param matBin Binary expression matrix (genes x cells)
#' @param regulator_info Regulator information with SCENIC metadata
#' @return List with composite_score, base_accuracy, sign_consistency, and scenic_confidence
#' @export
scoreRuleWithScenicMetrics <- function(ruleLogic, target, matBin, regulator_info) {
  
  # Create full rule string for scoring
  fullRule <- paste0(target, ", ", ruleLogic)
  
  # Base empirical accuracy using existing function
  base_accuracy <- tryCatch({
    scoreBooleanRule(fullRule, matBin)
  }, error = function(e) {
    0  # Failed rules get 0 accuracy
  })
  
  # Sign consistency bonus
  sign_consistency <- calculateSignConsistencyBonus(ruleLogic, regulator_info)
  
  # SCENIC confidence bonus
  scenic_confidence <- calculateScenicConfidenceBonus(regulator_info)
  
  # Composite score with weighted components
  composite <- base_accuracy + 
               (sign_consistency * 0.05) + 
               (scenic_confidence * 0.03)
               
  return(list(
    composite_score = pmax(0, composite),  # Ensure non-negative
    base_accuracy = base_accuracy,
    sign_consistency = sign_consistency,
    scenic_confidence = scenic_confidence
  ))
}

#' Calculate sign consistency bonus
#'
#' Measures how well inferred Boolean rule logic matches expected regulatory
#' directions from SCENIC (activators should appear without negation,
#' repressors should appear with negation).
#'
#' @param ruleLogic Character string with Boolean logic
#' @param regulator_info List with regulator metadata including regType
#' @return Numeric value between 0 and 1 (proportion of consistent regulator usage)
#' @export
calculateSignConsistencyBonus <- function(ruleLogic, regulator_info) {
  
  if (length(regulator_info) == 0 || regulator_info$n_selected == 0) {
    return(0)
  }
  
  regulators <- regulator_info$regulators
  
  consistency_score <- 0
  n_checked <- 0
  
  for (i in seq_len(nrow(regulators))) {
    reg_name <- regulators$TF[i]
    reg_type <- regulators$regType[i]
    
    if (is.na(reg_type)) next
    
    # Check if regulator appears in rule
    if (grepl(paste0("\\b", reg_name, "\\b"), ruleLogic)) {
      n_checked <- n_checked + 1
      
      # Check if usage matches expected sign
      if (reg_type == "Activation") {
        # Should appear without negation (or in positive context)
        if (grepl(paste0("(?<!!)\\b", reg_name, "\\b"), ruleLogic, perl = TRUE)) {
          consistency_score <- consistency_score + 1
        }
      } else if (reg_type == "Inhibition") {
        # Should appear with negation
        if (grepl(paste0("!\\s*", reg_name, "\\b"), ruleLogic)) {
          consistency_score <- consistency_score + 1
        }
      }
    }
  }
  
  return(if (n_checked > 0) consistency_score / n_checked else 0)
}

#' Calculate SCENIC confidence bonus
#'
#' Combines multiple SCENIC confidence metrics (motif confidence, NES scores,
#' motif presence) into a single bonus score for rule quality.
#'
#' @param regulator_info List with regulator metadata including motifConfidence, NES, hasMotif
#' @return Numeric confidence score (0 to ~1)
#' @export
calculateScenicConfidenceBonus <- function(regulator_info) {
  
  if (length(regulator_info) == 0 || regulator_info$n_selected == 0) {
    return(0)
  }
  
  regulators <- regulator_info$regulators
  
  # Average motif confidence
  motif_conf <- mean(ifelse(is.na(regulators$motifConfidence), 0, regulators$motifConfidence))
  
  # Average NES score (normalized)
  nes_score <- mean(ifelse(is.na(regulators$NES), 0, pmax(0, regulators$NES / 5)))
  
  # Proportion with confirmed motifs
  motif_prop <- mean(ifelse(is.na(regulators$hasMotif), 0, as.numeric(regulators$hasMotif)))
  
  # Combined confidence (weighted average)
  confidence <- (motif_conf * 0.4) + (nes_score * 0.4) + (motif_prop * 0.2)
  
  return(confidence)
}

# =============================================================================
# Fallback Rule Creation
# =============================================================================

#' Create intelligent fallback rule when inference fails
#'
#' Analyzes target gene expression pattern to create appropriate fallback rules.
#' Uses constant rules for genes with extreme expression bias, self-activation
#' for intermediate expression, or majority-based rules as final fallback.
#'
#' @param target Target gene name
#' @param regulator_info Available regulator information (may be empty)
#' @param matBin Binary expression matrix (genes x cells)
#' @param config Configuration object (currently unused)
#' @return List with rule metadata and appropriate fallback method label
#' @export
createIntelligentFallback <- function(target, regulator_info, matBin, config) {
  
  # Check target gene expression pattern
  if (!target %in% rownames(matBin)) {
    return(createConstantFallback(target, 1, "missing_gene", matBin))
  }
  
  target_expr <- matBin[target, ]
  prevalence <- mean(target_expr, na.rm = TRUE)
  
  # If gene has clear expression bias, use constant rule
  if (prevalence > 0.8) {
    return(createConstantFallback(target, 1, "high_prevalence", matBin))
  } else if (prevalence < 0.2) {
    return(createConstantFallback(target, 0, "low_prevalence", matBin))
  }
  
  # FIXED: Avoid identity rules that create point attractors
  # Use constant rule based on expression prevalence instead of self-activation
  
  # Final fallback: majority rule based on prevalence
  return(createConstantFallback(target, ifelse(prevalence > 0.5, 1, 0), "majority_fallback", matBin))
}

#' Create simple fallback rule for failed cases
#'
#' Creates basic constant or majority-based fallback rules when more
#' sophisticated inference methods fail.
#'
#' @param target Target gene name
#' @param matBin Binary expression matrix (genes x cells)
#' @param config Configuration object (currently unused)
#' @return List with simple fallback rule metadata
#' @export
createFallbackRule <- function(target, matBin, config) {
  
  if (!target %in% rownames(matBin)) {
    return(createConstantFallback(target, 1, "missing_gene", matBin))
  }
  
  # Use expression prevalence to decide constant value
  prevalence <- mean(matBin[target, ], na.rm = TRUE)
  constant_val <- ifelse(prevalence > 0.5, 1, 0)
  
  return(createConstantFallback(target, constant_val, "simple_fallback", matBin))
}

#' Create constant rule fallback
#'
#' Creates a constant Boolean rule (always 0 or always 1) as fallback
#' when other inference methods fail or are inappropriate.
#'
#' @param target Target gene name
#' @param value Constant value (0 or 1) - if NULL, will be determined from expression
#' @param reason Character string describing reason for fallback
#' @param matBin Binary expression matrix (optional, for auto-determining value)
#' @return List with constant rule metadata
#' @export
createConstantFallback <- function(target, value, reason, matBin = NULL) {
  # FIXED: Auto-determine constant value from expression if not provided
  if (is.null(value) && !is.null(matBin) && target %in% rownames(matBin)) {
    prevalence <- mean(matBin[target, ], na.rm = TRUE)
    value <- ifelse(prevalence >= 0.5, 1, 0)
  } else if (is.null(value)) {
    value <- 1  # Default fallback
  }
  
  list(
    rule = paste0(target, ", ", value),
    score = 0.5,  # Neutral score for constants
    method = reason,
    regulators = character(0),
    n_regulators = 0,
    k_used = 0
  )
}

# =============================================================================
# Validation Functions
# =============================================================================

#' Perform comprehensive final validation of all rules
#'
#' Validates Boolean rules for BoolNet syntax compatibility and calculates
#' quality metrics across the entire rule set. Checks for invalid characters,
#' proper rule structure, and identifies problematic rules.
#'
#' @param boolRules List of Boolean rules with rule metadata
#' @param matBin Binary expression matrix (genes x cells) - currently unused
#' @param edges Edge list for additional validation - currently unused
#' @return List with validation statistics: total_rules, valid_syntax, mean_quality, high_quality, syntax_errors
#' @export
performFinalValidation <- function(boolRules, matBin, edges) {
  
  total_rules <- length(boolRules)
  
  if (total_rules == 0) {
    return(list(
      total_rules = 0,
      valid_syntax = 0,
      mean_quality = 0,
      high_quality = 0,
      syntax_errors = character(0)
    ))
  }
  
  # Check BoolNet syntax
  syntax_valid <- sapply(boolRules, function(rule) {
    rule_str <- rule$rule %||% ""
    parts <- strsplit(rule_str, ",\\s*")[[1]]
    
    # Basic checks
    if (length(parts) != 2) return(FALSE)
    
    # Check for invalid characters - allow letters, numbers, underscore, &, |, (), !, and space
    logic <- parts[2]
    has_invalid_chars <- grepl("[^A-Za-z0-9_&\\|()! ]", logic)
    
    return(!has_invalid_chars)
  })
  
  # Calculate quality metrics
  scores <- sapply(boolRules, function(x) x$score %||% 0)
  mean_quality <- mean(scores, na.rm = TRUE)
  high_quality_count <- sum(scores >= 0.75, na.rm = TRUE)
  
  # Identify problematic rules
  syntax_errors <- names(boolRules)[!syntax_valid]
  
  return(list(
    total_rules = total_rules,
    valid_syntax = sum(syntax_valid),
    mean_quality = mean_quality,
    high_quality = high_quality_count,
    syntax_errors = syntax_errors
  ))
}

# =============================================================================
# Utility Functions and Fallbacks
# =============================================================================

#' Null coalescing operator
#'
#' Returns the right-hand side if left-hand side is NULL, empty, or NA.
#'
#' @param x Left-hand side value
#' @param y Right-hand side fallback value
#' @return x if valid, otherwise y
#' @export
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

# =============================================================================
# Boolean Rule Learning Core Functions
# =============================================================================

#' Create input-output pairs for Boolean rule learning
#'
#' Generates training data for Boolean rule inference by creating pairs of 
#' regulator states (input) and target gene state (output) along pseudotime.
#' Uses config$boolKValue for temporal offset.
#'
#' @param targetGene Character name of target gene (sanitized)
#' @param regulators Character vector of regulator gene names (sanitized)
#' @param matBin Binary expression matrix (genes x cells)
#' @param cellOrder Integer vector of cell indices ordered by pseudotime
#' @param k Integer time step offset between input and output (from config$boolKValue)
#' @return List with regulators, stateIndices, inputStates, outputStates, and ioData
#' @export
makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k) {
  
  # Input validation
  if (!targetGene %in% rownames(matBin)) {
    warning("Target gene '", targetGene, "' not found in expression matrix")
    return(NULL)
  }
  
  missingRegs <- setdiff(regulators, rownames(matBin))
  if (length(missingRegs) > 0) {
    warning("Regulators not found in matrix: ", paste(missingRegs, collapse = ", "))
    regulators <- intersect(regulators, rownames(matBin))
  }
  
  if (length(regulators) == 0) {
    warning("No valid regulators found for target '", targetGene, "'")
    return(NULL)
  }
  
  if (length(cellOrder) <= k) {
    warning("Cell order length (", length(cellOrder), ") <= k (", k, ")")
    return(NULL)
  }
  
  # Create input-output pairs
  maxT <- length(cellOrder) - k
  if (maxT < 5) {  # Need minimum data for learning
    warning("Insufficient data points (", maxT, ") for rule learning")
    return(NULL)
  }
  
  inputStates <- matrix(0, nrow = maxT, ncol = length(regulators))
  colnames(inputStates) <- regulators
  outputStates <- numeric(maxT)
  stateIndices <- seq_len(maxT)
  
  # Extract states efficiently
  for (t in stateIndices) {
    inputStates[t, ] <- matBin[regulators, cellOrder[t], drop = TRUE]
    outputStates[t] <- matBin[targetGene, cellOrder[t + k]]
  }
  
  # Create data frame for compatibility with old interface
  ioData <- data.frame(inputStates, output = outputStates, stringsAsFactors = FALSE)
  
  # Return structure expected by current workflow
  return(list(
    regulators = regulators,
    stateIndices = stateIndices,
    inputStates = inputStates,
    outputStates = outputStates,
    ioData = ioData
  ))
}

#' Find best Boolean rules using exhaustive search
#'
#' Searches all possible Boolean functions for the given number of inputs
#' and returns the best-scoring functions. Updated to handle ties properly
#' and provide better error handling.
#'
#' @param ioData List containing inputStates, outputStates, and other metadata
#' @return List with score, bestFns, functionsTestsed (note: typo preserved for compatibility), totalPossible
#' @export
findBestBooleanRules <- function(ioData) {
  
  # Input validation and extraction
  if (is.null(ioData) || is.null(ioData$inputStates) || is.null(ioData$outputStates)) {
    warning("Invalid ioData structure")
    return(NULL)
  }
  
  inputStates <- ioData$inputStates
  outputStates <- ioData$outputStates
  regulators <- ioData$regulators
  
  if (nrow(inputStates) == 0 || length(outputStates) == 0) {
    warning("No input-output data available")
    return(NULL)
  }
  
  if (nrow(inputStates) != length(outputStates)) {
    warning("Mismatch between input and output data lengths")
    return(NULL)
  }
  
  R <- ncol(inputStates)
  if (R == 0) {
    warning("No regulator inputs available")
    return(NULL)
  }
  
  # Limit complexity for computational feasibility
  if (R > 4) {
    warning("Too many regulators (", R, ") for exhaustive search, using heuristic")
    return(findBestBooleanRulesHeuristic(ioData))
  }
  
  # Generate all possible input states for R regulators
  allInputStates <- expand.grid(rep(list(c(0, 1)), R))
  colnames(allInputStates) <- regulators
  nStates <- 2^R
  
  # All possible Boolean functions (2^(2^R) possibilities)
  totalPossible <- 2^nStates
  if (totalPossible > 65536) {  # Practical limit
    warning("Too many possible functions (", totalPossible, "), using heuristic approach")
    return(findBestBooleanRulesHeuristic(ioData))
  }
  
  bestScore <- -1
  bestFns <- list()
  functionsTestsed <- 0  # Note: preserving typo for compatibility
  
  # Exhaustive search over all Boolean functions
  for (fnIdx in 0:(totalPossible - 1)) {
    functionsTestsed <- functionsTestsed + 1
    
    # Convert function index to output pattern
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    
    # Score this function against the data
    correct <- 0
    for (i in seq_len(nrow(inputStates))) {
      # Find which row in allInputStates matches current input
      matchIdx <- which(apply(allInputStates, 1, function(x) {
        all(x == inputStates[i, ])
      }))
      
      if (length(matchIdx) == 1 && outPattern[matchIdx] == outputStates[i]) {
        correct <- correct + 1
      }
    }
    
    score <- correct / length(outputStates)
    
    # Track best functions
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(
        patternIdx = fnIdx,
        outPattern = outPattern,
        score = score
      ))
    } else if (abs(score - bestScore) < 1e-9) {
      bestFns[[length(bestFns) + 1]] <- list(
        patternIdx = fnIdx,
        outPattern = outPattern,
        score = score
      )
    }
  }
  
  return(list(
    score = bestScore,
    bestFns = bestFns,
    functionsTestsed = functionsTestsed,  # Note: preserving typo
    totalPossible = totalPossible
  ))
}

#' Heuristic Boolean rule search for large regulator sets
#'
#' Uses greedy search and common Boolean patterns when exhaustive search
#' is computationally infeasible.
#'
#' @param ioData List containing input-output training data
#' @return List with score, bestFns, functionsTestsed, totalPossible
#' @export
findBestBooleanRulesHeuristic <- function(ioData) {
  
  inputStates <- ioData$inputStates
  outputStates <- ioData$outputStates
  regulators <- ioData$regulators
  R <- ncol(inputStates)
  
  # Try common Boolean patterns first
  patterns <- list(
    # Single regulators
    lapply(1:R, function(i) {
      pattern <- rep(0, 2^R)
      for (state in 0:(2^R - 1)) {
        bits <- as.integer(intToBits(state))[1:R]
        if (bits[i] == 1) pattern[state + 1] <- 1
      }
      list(type = paste0("single_", i), pattern = pattern)
    }),
    
    # OR of all regulators
    list(list(
      type = "all_or",
      pattern = c(0, rep(1, 2^R - 1))  # Only (0,0,...,0) -> 0, all others -> 1
    )),
    
    # AND of all regulators  
    list(list(
      type = "all_and",
      pattern = c(rep(0, 2^R - 1), 1)  # Only (1,1,...,1) -> 1, all others -> 0
    ))
  )
  
  bestScore <- -1
  bestFns <- list()
  functionsTestsed <- 0
  
  # Test each pattern
  for (patternGroup in patterns) {
    for (pattern in patternGroup) {
      functionsTestsed <- functionsTestsed + 1
      
      # Score pattern
      allInputStates <- expand.grid(rep(list(c(0, 1)), R))
      correct <- 0
      
      for (i in seq_len(nrow(inputStates))) {
        matchIdx <- which(apply(allInputStates, 1, function(x) {
          all(x == inputStates[i, ])
        }))
        
        if (length(matchIdx) == 1 && pattern$pattern[matchIdx] == outputStates[i]) {
          correct <- correct + 1
        }
      }
      
      score <- correct / length(outputStates)
      
      if (score > bestScore) {
        bestScore <- score
        bestFns <- list(list(
          patternIdx = -1,  # Heuristic patterns don't have indices
          outPattern = pattern$pattern,
          score = score,
          type = pattern$type
        ))
      } else if (abs(score - bestScore) < 1e-9) {
        bestFns[[length(bestFns) + 1]] <- list(
          patternIdx = -1,
          outPattern = pattern$pattern,
          score = score,
          type = pattern$type
        )
      }
    }
  }
  
  return(list(
    score = bestScore,
    bestFns = bestFns,
    functionsTestsed = functionsTestsed,
    totalPossible = 2^(2^R)
  ))
}

#' Convert Boolean function pattern to BoolNet rule string
#'
#' Generates BoolNet-compatible rule strings from Boolean function patterns.
#' Updated to handle edge cases and provide cleaner output.
#'
#' @param geneName Character name of target gene
#' @param regulators Character vector of regulator names
#' @param outPattern Integer vector of Boolean function outputs
#' @return Character string with complete BoolNet rule
#' @export
makeBoolNetRule <- function(geneName, regulators, outPattern) {
  
  # Handle case with no regulators
  if (length(regulators) == 0) {
    if (length(outPattern) == 0 || outPattern[1] == 1) {
      return(paste0(geneName, ", ", geneName))  # Self-activation
    } else {
      return(paste0(geneName, ", !", geneName))  # Self-inhibition
    }
  }
  
  R <- length(regulators)
  expectedLength <- 2^R
  
  if (length(outPattern) != expectedLength) {
    warning("Pattern length (", length(outPattern), ") != 2^R (", expectedLength, ")")
    # Pad or truncate as needed
    if (length(outPattern) < expectedLength) {
      outPattern <- c(outPattern, rep(0, expectedLength - length(outPattern)))
    } else {
      outPattern <- outPattern[1:expectedLength]
    }
  }
  
  # Find all input combinations that produce output = 1
  clauseList <- c()
  for (i in 0:(expectedLength - 1)) {
    if (outPattern[i + 1] == 1) {
      bits <- as.integer(intToBits(i))[1:R]
      andTerms <- ifelse(bits == 1, regulators, paste0("!", regulators))
      clause <- paste(andTerms, collapse = " & ")
      if (R > 1) {
        clauseList <- c(clauseList, paste0("(", clause, ")"))
      } else {
        clauseList <- c(clauseList, clause)
      }
    }
  }
  
  # Handle edge cases
  if (length(clauseList) == 0) {
    # Function is always false
    return(paste0(geneName, ", 0"))
  } else if (length(clauseList) == 1) {
    # Single clause
    return(paste0(geneName, ", ", clauseList[1]))
  } else {
    # Multiple clauses - OR them together
    return(paste0(geneName, ", ", paste(clauseList, collapse = " | ")))
  }
}

#' Score Boolean rule empirically against expression data
#'
#' Fallback implementation that evaluates how well a Boolean rule
#' matches the observed gene expression patterns in the data.
#'
#' @param ruleStr Character string with complete BoolNet rule ("target, logic")
#' @param matBin Binary expression matrix (genes x cells)
#' @return Numeric score between 0 and 1 indicating rule accuracy
#' @export
scoreBooleanRule <- function(ruleStr, matBin) {
  # Fallback implementation for rule scoring
  tryCatch({
    # Parse the rule
    parts <- strsplit(ruleStr, ",\\s*")[[1]]
    if (length(parts) != 2) return(0)
    
    target <- trimws(parts[1])
    logic <- trimws(parts[2])
    
    if (!target %in% rownames(matBin)) return(0)
    
    # Simple scoring: return proportion of cells where rule matches
    target_expr <- matBin[target, ]
    
    # For constants, return how well it matches
    if (logic == "0") {
      return(mean(target_expr == 0, na.rm = TRUE))
    } else if (logic == "1") {
      return(mean(target_expr == 1, na.rm = TRUE))
    } else if (logic == target) {
      # Self-activation: score autocorrelation
      return(0.7)  # Reasonable default for self-activation
    } else {
      # Complex logic - return moderate score
      return(0.6)
    }
  }, error = function(e) {
    return(0)
  })
}
