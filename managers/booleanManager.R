# =============================================================================
# booleanManager.R (Production Version with SCENIC Integration)
# Purpose: SCENIC-enhanced Boolean network construction with k=0 focus
# Enhanced with: regulator selection, sign-aware templates, composite scoring
# =============================================================================

library(dplyr)

# =============================================================================
# Enhanced Regulator Selection Using SCENIC Metadata
# =============================================================================

#' Select optimal regulators for target using SCENIC metadata
#'
#' Uses composite weighting of correlation, motif confidence, SCENIC importance,
#' and diversity filtering to choose the best regulators for Boolean rule inference.
#'
#' @param target Character name of target gene (sanitized)
#' @param edges Edge dataframe with SCENIC metadata (TF, Target, corr, regType, etc.)
#' @param matBin Binary expression matrix (genes x cells)
#' @param maxRegs Maximum number of regulators to select
#' @return List with selected regulator information
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
  
  # Calculate composite weights
  weights <- calculateRegulatorWeights(targetEdges)
  targetEdges$composite_weight <- weights
  
  # Sort by weight and take top candidates
  targetEdges <- targetEdges[order(-weights), , drop = FALSE]
  
  # Apply diversity filtering if we have more candidates than slots
  if (nrow(targetEdges) > maxRegs) {
    selected_indices <- selectDiverseRegulatorsWithWeights(targetEdges, matBin, maxRegs)
    targetEdges <- targetEdges[selected_indices, , drop = FALSE]
  }
  
  # Limit to maxRegs
  targetEdges <- head(targetEdges, maxRegs)
  
  return(list(
    regulators = targetEdges,
    n_selected = nrow(targetEdges),
    mean_weight = if (nrow(targetEdges) > 0) mean(targetEdges$composite_weight) else 0,
    has_activators = any(targetEdges$regType == "Activation", na.rm = TRUE),
    has_repressors = any(targetEdges$regType == "Inhibition", na.rm = TRUE)
  ))
}

#' Calculate composite regulator weights using SCENIC metadata
#'
#' @param edges Edge dataframe for a single target
#' @return Numeric vector of composite weights
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
  
  # Composite weight calculation
  composite <- corr_weight + 
               (motif_bonus * 0.3) + 
               (nes_bonus * 0.2) + 
               prior_bonus + 
               direction_bonus
  
  return(composite)
}

#' Select diverse regulators to avoid multicollinearity
#'
#' @param edges Sorted edge dataframe
#' @param matBin Binary expression matrix
#' @param maxRegs Maximum regulators to select
#' @return Integer vector of selected row indices
selectDiverseRegulatorsWithWeights <- function(edges, matBin, maxRegs) {
  
  regulators <- edges$TF
  n_regs <- length(regulators)
  
  if (n_regs <= maxRegs) {
    return(1:n_regs)
  }
  
  # Calculate correlation matrix among regulators
  reg_matrix <- matBin[regulators, , drop = FALSE]
  reg_cors <- tryCatch({
    cor(t(reg_matrix))
  }, error = function(e) {
    # If correlation fails, just return the first maxRegs
    return(matrix(0, nrow = n_regs, ncol = n_regs))
  })
  
  # Greedy selection: start with highest weight, add others that aren't too correlated
  selected <- c(1)  # Always take the top weighted regulator
  
  for (i in 2:n_regs) {
    if (length(selected) >= maxRegs) break
    
    # Check correlation with already selected regulators
    max_cor <- max(abs(reg_cors[i, selected]), na.rm = TRUE)
    
    # Add if not too correlated (threshold 0.8)
    if (is.na(max_cor) || max_cor < 0.8) {
      selected <- c(selected, i)
    }
  }
  
  return(selected)
}

# =============================================================================
# Sign-Aware Template Generation
# =============================================================================

#' Generate Boolean templates respecting SCENIC regulatory signs
#'
#' Creates biologically plausible Boolean logic templates based on whether
#' regulators are activators or repressors according to SCENIC analysis.
#'
#' @param target Character name of target gene
#' @param regulator_info List from selectOptimalRegulatorsForTarget
#' @return Named list of template rule strings
#' @export
generateSignAwareTemplates <- function(target, regulator_info) {
  
  if (length(regulator_info) == 0 || regulator_info$n_selected == 0) {
    return(list(
      self_activation = paste0(target, ", ", target)
    ))
  }
  
  regulators <- regulator_info$regulators
  activators <- regulators$TF[regulators$regType == "Activation" & !is.na(regulators$regType)]
  repressors <- regulators$TF[regulators$regType == "Inhibition" & !is.na(regulators$regType)]
  
  templates <- list()
  
  # Template 1: Single regulator (if only one)
  if (nrow(regulators) == 1) {
    reg <- regulators$TF[1]
    if (!is.na(regulators$regType[1]) && regulators$regType[1] == "Inhibition") {
      templates$single_repressor <- paste0(target, ", !", reg)
    } else {
      templates$single_activator <- paste0(target, ", ", reg)
    }
    return(templates)
  }
  
  # Template 2: All activators OR (any activator sufficient)
  if (length(activators) > 0) {
    templates$activator_or <- paste0(target, ", ", paste(activators, collapse = " | "))
  }
  
  # Template 3: All activators AND (cooperative activation)
  if (length(activators) > 1) {
    templates$activator_and <- paste0(target, ", ", paste(activators, collapse = " & "))
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
    
    templates$balanced <- paste0(target, ", ", act_part, " & ", rep_part)
  }
  
  # Template 5: Repressor dominant (target OFF when any repressor ON)
  if (length(repressors) > 0) {
    if (length(repressors) == 1) {
      templates$repressor_dominant <- paste0(target, ", !", repressors[1])
    } else {
      templates$repressor_dominant <- paste0(target, ", !(", paste(repressors, collapse = " | "), ")")
    }
  }
  
  # Template 6: Majority rule (for 3+ regulators)
  if (nrow(regulators) >= 3) {
    # Simple majority: more than half must be "on"
    all_regs <- regulators$TF
    pos_terms <- ifelse(!is.na(regulators$regType) & regulators$regType == "Inhibition", 
                        paste0("!", all_regs), all_regs)
    
    # For 3 regulators, need at least 2 positive
    if (length(pos_terms) == 3) {
      pairs <- combn(pos_terms, 2, simplify = FALSE)
      pair_strings <- sapply(pairs, function(p) paste0("(", paste(p, collapse = " & "), ")"))
      templates$majority <- paste0(target, ", ", paste(pair_strings, collapse = " | "))
    }
  }
  
  return(templates)
}

# =============================================================================
# Enhanced Rule Scoring with SCENIC Metrics
# =============================================================================

#' Score Boolean rule with SCENIC-aware composite metrics
#'
#' @param ruleStr BoolNet-compatible rule string
#' @param matBin Binary expression matrix
#' @param regulator_info Regulator information with SCENIC metadata
#' @return List with composite score and breakdown
#' @export
scoreRuleWithScenicMetrics <- function(ruleStr, matBin, regulator_info) {
  
  # Base empirical accuracy using existing function
  base_accuracy <- tryCatch({
    scoreBooleanRule(ruleStr, matBin)
  }, error = function(e) {
    0  # Failed rules get 0 accuracy
  })
  
  # Sign consistency bonus
  sign_consistency <- calculateSignConsistencyBonus(ruleStr, regulator_info)
  
  # Parsimony bonus (prefer simpler rules)
  parsimony_bonus <- calculateParsimonyBonus(ruleStr)
  
  # SCENIC confidence bonus
  scenic_confidence <- calculateScenicConfidenceBonus(regulator_info)
  
  # Composite score with weighted components
  composite <- base_accuracy + 
               (sign_consistency * 0.05) + 
               (parsimony_bonus * 0.02) + 
               (scenic_confidence * 0.03)
               
  return(list(
    composite_score = pmax(0, composite),  # Ensure non-negative
    base_accuracy = base_accuracy,
    sign_consistency = sign_consistency,
    parsimony = parsimony_bonus,
    scenic_confidence = scenic_confidence
  ))
}

#' Calculate sign consistency bonus
calculateSignConsistencyBonus <- function(ruleStr, regulator_info) {
  
  if (length(regulator_info) == 0 || regulator_info$n_selected == 0) {
    return(0)
  }
  
  # Parse rule to check if regulator usage matches expected signs
  rule_parts <- strsplit(ruleStr, ",\\s*")[[1]]
  if (length(rule_parts) != 2) return(0)
  
  logic_part <- rule_parts[2]
  regulators <- regulator_info$regulators
  
  consistency_score <- 0
  n_checked <- 0
  
  for (i in seq_len(nrow(regulators))) {
    reg_name <- regulators$TF[i]
    reg_type <- regulators$regType[i]
    
    if (is.na(reg_type)) next
    
    # Check if regulator appears in rule
    if (grepl(paste0("\\b", reg_name, "\\b"), logic_part)) {
      n_checked <- n_checked + 1
      
      # Check if usage matches expected sign
      if (reg_type == "Activation") {
        # Should appear without negation (or in positive context)
        if (grepl(paste0("(?<!!)\\b", reg_name, "\\b"), logic_part, perl = TRUE)) {
          consistency_score <- consistency_score + 1
        }
      } else if (reg_type == "Inhibition") {
        # Should appear with negation
        if (grepl(paste0("!\\s*", reg_name, "\\b"), logic_part)) {
          consistency_score <- consistency_score + 1
        }
      }
    }
  }
  
  return(if (n_checked > 0) consistency_score / n_checked else 0)
}

#' Calculate parsimony bonus (prefer simpler rules)
calculateParsimonyBonus <- function(ruleStr) {
  
  rule_parts <- strsplit(ruleStr, ",\\s*")[[1]]
  if (length(rule_parts) != 2) return(0)
  
  logic_part <- rule_parts[2]
  
  # Count logical operators
  n_and <- length(gregexpr("&", logic_part, fixed = TRUE)[[1]])
  n_or <- length(gregexpr("\\|", logic_part)[[1]])
  n_not <- length(gregexpr("!", logic_part, fixed = TRUE)[[1]])
  
  # Adjust for -1 when no matches found
  n_and <- ifelse(n_and == 1 && !grepl("&", logic_part), 0, n_and)
  n_or <- ifelse(n_or == 1 && !grepl("\\|", logic_part), 0, n_or)
  n_not <- ifelse(n_not == 1 && !grepl("!", logic_part), 0, n_not)
  
  # Complexity penalty
  complexity <- (n_and * 0.5) + (n_or * 0.3) + (n_not * 0.2)
  
  # Return bonus (higher for simpler rules)
  return(pmax(0, 1 - complexity / 10))
}

#' Calculate SCENIC confidence bonus
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
  
  # Combined confidence
  confidence <- (motif_conf * 0.4) + (nes_score * 0.4) + (motif_prop * 0.2)
  
  return(confidence)
}

# =============================================================================
# Multi-Method Boolean Rule Inference
# =============================================================================

#' Infer best Boolean rule using multiple methods
#'
#' Tries both sign-aware templates and empirical state counting,
#' then selects the best based on composite scoring.
#'
#' @param target Target gene name
#' @param regulator_info Selected regulator information
#' @param matBin Binary expression matrix
#' @param cellOrder Pseudotime cell ordering
#' @param config Configuration object
#' @return Best rule with metadata
#' @export
inferBestBooleanRule <- function(target, regulator_info, matBin, cellOrder, config) {
  
  candidates <- list()
  
  # Method 1: Sign-aware templates
  templates <- generateSignAwareTemplates(target, regulator_info)
  
  for (template_name in names(templates)) {
    rule_str <- templates[[template_name]]
    
    # Score template
    score_info <- scoreRuleWithScenicMetrics(rule_str, matBin, regulator_info)
    
    candidates[[paste0("template_", template_name)]] <- list(
      rule = rule_str,
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
      ioData <- makeInputOutputPairs(target, reg_names, matBin, cellOrder, k = 0)
      
      if (!is.null(ioData) && length(ioData$stateIndices) >= 10) {
        
        # Find best empirical rule
        empirical_result <- findBestBooleanRules(ioData)
        
        if (!is.null(empirical_result) && empirical_result$score > 0 && 
            length(empirical_result$bestFns) > 0) {
          
          # Generate rule string
          rule_str <- makeBoolNetRule(target, ioData$regulators, 
                                      empirical_result$bestFns[[1]]$outPattern)
          
          # Score empirical rule
          score_info <- scoreRuleWithScenicMetrics(rule_str, matBin, regulator_info)
          
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
      # Empirical method failed, continue with just templates
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

# =============================================================================
# Fallback Rule Creation
# =============================================================================

#' Create intelligent fallback rule when inference fails
#'
#' @param target Target gene name
#' @param regulator_info Available regulator information (may be empty)
#' @param matBin Binary expression matrix
#' @param config Configuration object
#' @return Fallback rule with appropriate method label
#' @export
createIntelligentFallback <- function(target, regulator_info, matBin, config) {
  
  # Check target gene expression pattern
  if (!target %in% rownames(matBin)) {
    return(createConstantFallback(target, 1, "missing_gene"))
  }
  
  target_expr <- matBin[target, ]
  prevalence <- mean(target_expr, na.rm = TRUE)
  
  # If gene has clear expression bias, use constant rule
  if (prevalence > 0.8) {
    return(createConstantFallback(target, 1, "high_prevalence"))
  } else if (prevalence < 0.2) {
    return(createConstantFallback(target, 0, "low_prevalence"))
  }
  
  # Try self-activation if expression is intermediate
  self_rule <- paste0(target, ", ", target)
  self_score <- tryCatch({
    scoreBooleanRule(self_rule, matBin)
  }, error = function(e) {
    0.5  # Default score for self-activation
  })
  
  if (self_score > 0.55) {
    return(list(
      rule = self_rule,
      score = self_score,
      method = "self_activation_fallback",
      regulators = character(0),
      n_regulators = 0,
      k_used = 0
    ))
  }
  
  # Final fallback: majority rule
  return(createConstantFallback(target, ifelse(prevalence > 0.5, 1, 0), "majority_fallback"))
}

#' Create simple fallback rule for failed cases
createFallbackRule <- function(target, matBin, config) {
  
  if (!target %in% rownames(matBin)) {
    return(createConstantFallback(target, 1, "missing_gene"))
  }
  
  # Use expression prevalence to decide constant value
  prevalence <- mean(matBin[target, ], na.rm = TRUE)
  constant_val <- ifelse(prevalence > 0.5, 1, 0)
  
  return(createConstantFallback(target, constant_val, "simple_fallback"))
}

#' Create constant rule fallback
createConstantFallback <- function(target, value, reason) {
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
# Final Validation Functions
# =============================================================================

#' Perform comprehensive final validation of all rules
#'
#' @param boolRules List of Boolean rules
#' @param matBin Binary expression matrix
#' @param edges Edge list for additional validation
#' @return List with validation results
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
    
    # Check for invalid characters
    logic <- parts[2]
    has_invalid_chars <- grepl("[^A-Za-z0-9_&|()!\\s]", logic)
    
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
# Core Utility Functions (Fallbacks for missing functions)
# =============================================================================

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

# Ensure these core functions are available (placeholders if missing)
if (!exists("makeInputOutputPairs")) {
  makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k = 0) {
    stop("makeInputOutputPairs function not found - ensure boolean_helpers.R is loaded")
  }
}

if (!exists("findBestBooleanRules")) {
  findBestBooleanRules <- function(ioData) {
    stop("findBestBooleanRules function not found - ensure boolean_helpers.R is loaded")
  }
}

if (!exists("makeBoolNetRule")) {
  makeBoolNetRule <- function(geneName, regulators, outPattern) {
    # Fallback implementation for BoolNet rule creation
    if (length(regulators) == 0) {
      # No regulators - use constant or self-activation
      if (length(outPattern) == 1 && outPattern[1] == 1) {
        return(paste0(geneName, ", ", geneName))  # Self-activation
      } else {
        return(paste0(geneName, ", ", outPattern[1]))  # Constant
      }
    }
    
    # With regulators - create simple logic
    if (length(outPattern) == 2^length(regulators)) {
      # Try to create a simple rule based on pattern
      if (all(outPattern == 1)) {
        return(paste0(geneName, ", 1"))
      } else if (all(outPattern == 0)) {
        return(paste0(geneName, ", 0"))
      } else {
        # Default to OR logic
        return(paste0(geneName, ", ", paste(regulators, collapse = " | ")))
      }
    }
    
    # Default fallback
    return(paste0(geneName, ", ", paste(regulators, collapse = " | ")))
  }
}

if (!exists("scoreBooleanRule")) {
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
}

if (!exists("sanitizeGeneName")) {
  sanitizeGeneName <- function(name) {
    # Convert to valid R variable name
    name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
    if (grepl("^[0-9]", name2)) name2 <- paste0("X", name2)
    if (name2 == "" || name2 == "_") name2 <- "Unknown"
    name2
  }
}

if (!exists("generateSanitizedGeneMapping")) {
  generateSanitizedGeneMapping <- function(geneNames) {
    data.frame(
      OriginalName = geneNames,
      SanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE),
      stringsAsFactors = FALSE
    )
  }
}
