# =============================================================================
# booleanManager.R
# Purpose: Boolean network construction, rule inference, and optimization
# Contains both template-based and exhaustive Boolean rule inference methods
# =============================================================================

# =============================================================================
# Template-Based Boolean Rule Synthesis
# =============================================================================

#' Generate multiple Boolean rule templates for a target gene
#'
#' Creates various Boolean logic templates (cooperative, redundant, majority, etc.)
#' for a target gene based on its regulators and their regulatory signs.
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector (1 for activation, -1 for repression)
#' @return List of rule templates with descriptions
#' @export
generateBooleanTemplates <- function(targetGene, regulators, regulatorSigns) {
  
  if (length(regulators) == 0) {
    return(list(
      self_activation = list(
        rule = paste0(targetGene, ", ", targetGene),
        description = "Self-activation (default state)"
      )
    ))
  }
  
  # Prepare terms with signs
  activators <- regulators[regulatorSigns[regulators] > 0]
  repressors <- regulators[regulatorSigns[regulators] < 0]
  
  templates <- list()
  
  # 1. SINGLE REGULATOR TEMPLATES (if only 1 regulator)
  if (length(regulators) == 1) {
    reg <- regulators[1]
    if (regulatorSigns[reg] > 0) {
      templates$single_activator <- list(
        rule = paste0(targetGene, ", ", reg),
        description = paste("Single activator:", reg)
      )
    } else {
      templates$single_repressor <- list(
        rule = paste0(targetGene, ", !", reg),
        description = paste("Single repressor:", reg)
      )
    }
    return(templates)
  }
  
  # 2. COOPERATIVE ACTIVATION (AND logic)
  if (length(activators) >= 2) {
    andTerms <- c(activators, if(length(repressors) > 0) paste0("!", repressors) else character(0))
    templates$cooperative_activation <- list(
      rule = paste0(targetGene, ", ", paste(andTerms, collapse = " & ")),
      description = "Cooperative activation (all activators + no repressors)"
    )
  }
  
  # 3. REDUNDANT ACTIVATION (OR logic for activators)
  if (length(activators) >= 2) {
    # Any activator works, but repressors can block
    if (length(repressors) > 0) {
      orActivators <- paste0("(", paste(activators, collapse = " | "), ")")
      andRepressors <- paste0("!(", paste(repressors, collapse = " | "), ")")
      ruleLogic <- paste(orActivators, andRepressors, sep = " & ")
    } else {
      ruleLogic <- paste(activators, collapse = " | ")
    }
    
    templates$redundant_activation <- list(
      rule = paste0(targetGene, ", ", ruleLogic),
      description = "Redundant activation (any activator works)"
    )
  }
  
  # 4. BALANCED ACTIVATION (mix of AND and OR)
  if (length(activators) >= 2 && length(repressors) >= 1) {
    # Need at least one activator AND no repressors
    orActivators <- paste0("(", paste(activators, collapse = " | "), ")")
    andRepressors <- paste0("!", paste(repressors, collapse = " & !"))
    
    templates$balanced_regulation <- list(
      rule = paste0(targetGene, ", ", orActivators, " & ", andRepressors),
      description = "Balanced: any activator + all repressors off"
    )
  }
  
  # 5. MAJORITY RULE (for 3+ regulators)
  if (length(regulators) >= 3) {
    # Simple majority: more activators than repressors
    majorityTerms <- character()
    
    # All possible combinations where activators outnumber active repressors
    if (length(activators) >= 2) {
      # At least 2 activators on
      activatorPairs <- combn(activators, 2, simplify = FALSE)
      for (pair in activatorPairs) {
        pairTerm <- paste(pair, collapse = " & ")
        majorityTerms <- c(majorityTerms, pairTerm)
      }
      
      # If we have 3+ activators, add triple combinations
      if (length(activators) >= 3) {
        allActivators <- paste(activators, collapse = " & ")
        majorityTerms <- c(majorityTerms, allActivators)
      }
    }
    
    if (length(majorityTerms) > 0) {
      templates$majority_rule <- list(
        rule = paste0(targetGene, ", ", paste(unique(majorityTerms), collapse = " | ")),
        description = "Majority rule (multiple activator combinations)"
      )
    }
  }
  
  # 6. SIMPLE AND (default fallback)
  allTerms <- c(activators, if(length(repressors) > 0) paste0("!", repressors) else character(0))
  templates$simple_and <- list(
    rule = paste0(targetGene, ", ", paste(allTerms, collapse = " & ")),
    description = "Simple AND (all activators + no repressors)"
  )
  
  # 7. DOMINANT ACTIVATOR (strongest activator alone)
  if (length(activators) >= 1) {
    # Assume first activator is strongest (could sort by SCENIC score)
    dominantActivator <- activators[1]
    templates$dominant_activator <- list(
      rule = paste0(targetGene, ", ", dominantActivator),
      description = paste("Dominant activator:", dominantActivator)
    )
  }
  
  return(templates)
}

#' Test multiple templates and select the best one
#'
#' Evaluates multiple Boolean rule templates against expression data and
#' returns the template with the highest accuracy score.
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector of regulatory signs
#' @param matBin Binary expression matrix (genes x cells)
#' @return List with best rule, score, and all tested templates
#' @export
synthesizeBestBooleanRule <- function(targetGene, regulators, regulatorSigns, matBin) {
  
  # Generate all possible templates
  templates <- generateBooleanTemplates(targetGene, regulators, regulatorSigns)
  
  # Score each template
  templateScores <- list()
  bestScore <- -1
  bestTemplate <- NULL
  bestRuleName <- NULL
  
  for (templateName in names(templates)) {
    template <- templates[[templateName]]
    ruleStr <- template$rule
    
    tryCatch({
      score <- scoreBooleanRule(ruleStr, matBin)
      templateScores[[templateName]] <- list(
        rule = ruleStr,
        score = score,
        description = template$description
      )
      
      if (score > bestScore) {
        bestScore <- score
        bestTemplate <- template
        bestRuleName <- templateName
      }
    }, error = function(e) {
      # Skip templates that fail to evaluate
      templateScores[[templateName]] <- list(
        rule = ruleStr,
        score = 0,
        description = paste(template$description, "(FAILED)")
      )
    })
  }
  
  return(list(
    bestRule = bestTemplate$rule,
    bestScore = bestScore,
    bestTemplateName = bestRuleName,
    bestDescription = bestTemplate$description,
    allTemplates = templateScores,
    regulators = regulators,
    nRegulators = length(regulators)
  ))
}

#' Enhanced batch synthesis with template variety
#'
#' Processes multiple genes and synthesizes Boolean rules using predefined
#' templates. Tests various logic patterns and selects the best-fitting template.
#'
#' @param edges Data frame with columns TF, Target, regType, corr
#' @param matBin Binary expression matrix (genes x cells)
#' @param maxRegulators Maximum number of regulators per target gene
#' @return Named list of rule objects with template information
#' @export
synthesizeBooleanRulesBatch <- function(edges, matBin, maxRegulators = 5) {
  
  targets <- unique(edges$Target)
  total <- length(targets)
  
  if (total == 0) {
    warning("No targets found in edges data frame")
    return(list())
  }
  
  rules <- list()
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in seq_along(targets)) {
    gene <- targets[i]
    setTxtProgressBar(pb, i)
    
    # Get edges for this gene
    subEdges <- edges[edges$Target == gene, ]
    if (nrow(subEdges) == 0) {
      next
    }
    
    # Get regulators
    topEdges <- head(subEdges[order(-abs(subEdges$corr)), ], maxRegulators)
    regulators <- topEdges$TF
    
    # Check if gene and regulators are in matrix
    geneInMatrix <- gene %in% rownames(matBin)
    regsInMatrix <- regulators %in% rownames(matBin)
    
    if (!geneInMatrix || !any(regsInMatrix)) {
      next
    }
    
    # Create signs
    signs <- setNames(ifelse(topEdges$regType == "Activation", 1, -1), topEdges$TF)
    
    # Try to synthesize rule
    tryCatch({
      ruleResult <- synthesizeBestBooleanRule(gene, regulators, signs, matBin)
      
      rules[[gene]] <- list(
        rule = ruleResult$bestRule,
        regulators = regulators,
        bestScore = ruleResult$bestScore,
        score = ruleResult$bestScore,
        nRegulators = length(regulators),
        templateUsed = ruleResult$bestTemplateName,
        templateDescription = ruleResult$bestDescription,
      )
      
    }, error = function(e) {
      # Silently skip genes that fail rule synthesis
      next
    })
  }
  
  close(pb)
  return(rules)
}

#' Helper function to score Boolean rules
#'
#' Evaluates a Boolean rule against binary expression data and returns
#' an accuracy score (proportion of correct predictions).
#'
#' @param ruleStr BoolNet-style rule string (e.g., "Gene1, Gene2 & !Gene3")
#' @param matBin Binary expression matrix (genes x cells)
#' @return Accuracy score between 0 and 1
#' @export
scoreBooleanRule <- function(ruleStr, matBin) {
  
  # Parse rule
  ruleParts <- strsplit(ruleStr, ",\\s*")[[1]]
  if (length(ruleParts) != 2) {
    warning("Invalid rule format: ", ruleStr)
    return(0)
  }
  
  targetGene <- ruleParts[1]
  logic <- ruleParts[2]
  
  if (!targetGene %in% rownames(matBin)) {
    warning("Target gene not found in matrix: ", targetGene)
    return(0)
  }
  
  n <- ncol(matBin)
  predictions <- logical(n)
  actual <- as.logical(matBin[targetGene, ])
  
  # Check if all required genes are in the matrix
  logicCleaned <- gsub("[&|()!]", " ", logic)
  requiredGenes <- unique(trimws(unlist(strsplit(logicCleaned, "\\s+"))))
  requiredGenes <- requiredGenes[requiredGenes != ""]
  requiredGenes <- requiredGenes[requiredGenes %in% rownames(matBin)]
  
  if (length(requiredGenes) == 0) {
    warning("No valid regulator genes found for rule: ", ruleStr)
    return(0)
  }
  
  # Evaluate rule for each cell
  for (i in 1:n) {
    env <- as.list(as.logical(matBin[, i]))
    names(env) <- rownames(matBin)
    
    predictions[i] <- tryCatch({
      result <- eval(parse(text = logic), envir = env)
      if (is.logical(result) && length(result) == 1 && !is.na(result)) {
        result
      } else {
        FALSE
      }
    }, error = function(e) {
      FALSE
    })
  }
  
  score <- mean(predictions == actual, na.rm = TRUE)
  return(score)
}

# =============================================================================
# Exhaustive Boolean Rule Learning (k-step approach)
# =============================================================================

#' Compute input-output pairs for Boolean rule learning (optimized)
#'
#' Creates input-output pairs by taking regulator states at time t and
#' target gene state at time t+k. Optimized version that avoids per-row
#' data frame construction and returns compact state indices.
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param matBin Binary expression matrix (genes x cells)
#' @param cellOrder Integer vector specifying cell ordering (e.g., by pseudotime)
#' @param k Integer step size for temporal relationship (0 = same time)
#' @return List with stateIndices, outputs, and regulators
#' @export
makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k = 0) {
  
  # Validate inputs
  if (!targetGene %in% rownames(matBin)) {
    warning("Target gene ", targetGene, " not found in matrix")
    return(NULL)
  }
  
  missingRegs <- regulators[!regulators %in% rownames(matBin)]
  if (length(missingRegs) > 0) {
    warning("Regulators not found in matrix: ", paste(missingRegs, collapse = ", "))
    regulators <- regulators[regulators %in% rownames(matBin)]
  }
  
  if (length(regulators) == 0) {
    warning("No valid regulators found for ", targetGene)
    return(NULL)
  }
  
  maxT <- length(cellOrder) - k
  
  if (maxT <= 0) {
    warning("Not enough cells for k = ", k, " with ", length(cellOrder), " ordered cells")
    return(NULL)
  }
  
  # Check for constant target gene (optimization)
  targetCells <- cellOrder[seq_len(maxT) + k]
  targetVals  <- matBin[targetGene, targetCells]
  if (length(unique(targetVals)) == 1) {
    warning("Target gene ", targetGene, " is constant - skipping")
    return(NULL)
  }
  
  # Extract regulator matrix for all relevant timepoints
  inputCells <- cellOrder[seq_len(maxT)]
  regMatrix  <- matBin[regulators, inputCells, drop = FALSE]
  
  # Remove constant regulators (optimization)
  regVariance <- apply(regMatrix, 1, function(x) length(unique(x)))
  variableRegs <- regulators[regVariance > 1]
  
  if (length(variableRegs) == 0) {
    warning("All regulators are constant for ", targetGene, " - skipping")
    return(NULL)
  }
  
  if (length(variableRegs) < length(regulators)) {
    #message("Removed ", length(regulators) - length(variableRegs), " constant regulators for ", targetGene)
    regulators <- variableRegs
    regMatrix <- regMatrix[regulators, , drop = FALSE]
  }
  
  R <- length(regulators)
  
  # Convert regulator states to state indices using binary encoding
  weights <- 2^(0:(R-1))
  stateIndices <- 1 + as.integer(regMatrix[1, ] * weights[1])
  
  if (R > 1) {
    for (i in 2:R) {
      stateIndices <- stateIndices + as.integer(regMatrix[i, ] * weights[i])
    }
  }
  
  # Get target outputs
  outputs <- as.integer(targetVals)
  
  return(list(
    stateIndices = stateIndices,
    outputs = outputs,
    regulators = regulators,
    nStates = 2^R
  ))
}

#' Find the best Boolean rules using optimized state counting
#'
#' Uses state counting instead of exhaustive enumeration to find optimal
#' Boolean functions. Mathematically equivalent to exhaustive search with
#' tie-breaking by OR, but orders of magnitude faster.
#'
#' @param ioData List from makeInputOutputPairs with stateIndices, outputs, etc.
#' @return List with best score and best functions
#' @export
findBestBooleanRules <- function(ioData) {
  
  if (is.null(ioData)) {
    warning("No input-output data provided")
    return(list(score = 0, bestFns = list()))
  }
  
  stateIndices <- ioData$stateIndices
  outputs <- ioData$outputs
  regulators <- ioData$regulators
  nStates <- ioData$nStates
  
  if (length(stateIndices) == 0 || length(outputs) == 0) {
    warning("Empty state indices or outputs")
    return(list(score = 0, bestFns = list()))
  }
  
  # Count outcomes for each unique state
  count1 <- integer(nStates)  # Count of output=1 for each state
  count0 <- integer(nStates)  # Count of output=0 for each state
  
  for (i in seq_along(stateIndices)) {
    stateIdx <- stateIndices[i]
    if (outputs[i] == 1) {
      count1[stateIdx] <- count1[stateIdx] + 1
    } else {
      count0[stateIdx] <- count0[stateIdx] + 1
    }
  }
  
  # For each state, choose the majority outcome
  # If tied, choose 1 (equivalent to OR-ing tied best functions)
  outPattern <- as.integer(count1 >= count0)
  
  # Calculate the score of this optimal function
  correct <- 0
  total <- length(stateIndices)
  
  for (i in seq_along(stateIndices)) {
    stateIdx <- stateIndices[i]
    if (outPattern[stateIdx] == outputs[i]) {
      correct <- correct + 1
    }
  }
  
  score <- correct / total
  
  # Create a function object in the same format as the old approach
  bestFn <- list(
    patternIdx = 0,  # Not meaningful in this approach
    outPattern = outPattern,
    score = score
  )
  
  return(list(
    score = score,
    bestFns = list(bestFn),
    functionsTested = nStates,  # We "tested" one function per state - fixed typo
    totalPossible = 2^nStates,  # What exhaustive would have tested
    earlyStop = FALSE,  # Not applicable to this method
    method = "state_counting"
  ))
}

#' Combine multiple Boolean functions using logical OR
#'
#' Takes a list of Boolean functions and combines them using OR logic.
#' Used when multiple functions tie for the best score.
#'
#' @param fnList List of Boolean functions from findBestBooleanRules
#' @return Integer vector representing combined Boolean function output
#' @export
combineBooleanFunctionsByOr <- function(fnList) {
  if (length(fnList) == 0) {
    warning("Empty function list provided")
    return(integer(0))
  }
  
  if (length(fnList) == 1) {
    return(fnList[[1]]$outPattern)
  }
  
  # Combine multiple functions with OR
  mat <- sapply(fnList, `[[`, "outPattern")
  return(as.integer(rowSums(mat) > 0))
}

# =============================================================================
# Rule String Generation and Utilities
# =============================================================================

#' Generate BoolNet-compatible rule string from Boolean function
#'
#' Converts a Boolean function output pattern into a human-readable rule string
#' suitable for BoolNet network specification.
#'
#' @param geneName Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param outPattern Integer vector representing Boolean function output
#' @return Character string in BoolNet rule format
#' @export
makeBoolNetRule <- function(geneName, regulators, outPattern) {
  
  # Handle case with no regulators
  if (length(regulators) == 0) {
    ruleStr <- ifelse(outPattern[1] == 1, geneName, paste0("!", geneName))
    return(paste0(geneName, ", ", ruleStr))
  }
  
  R <- length(regulators)
  expectedLength <- 2^R
  
  if (length(outPattern) != expectedLength) {
    stop("Output pattern length (", length(outPattern),
         ") doesn't match expected length (", expectedLength, ") for ", R, " regulators")
  }
  
  # Build DNF (Disjunctive Normal Form) from truth table
  clauseList <- character()
  
  for (i in 0:(expectedLength - 1)) {
    if (outPattern[i + 1] == 1) {
      # Convert index to binary representation
      bits <- as.integer(intToBits(i))[1:R]
      andTerms <- ifelse(bits == 1, regulators, paste0("!", regulators))
      clause <- paste0("(", paste(andTerms, collapse = " & "), ")")
      clauseList <- c(clauseList, clause)
    }
  }
  
  # Handle case where function is always false
  if (length(clauseList) == 0) {
    return(paste0(geneName, ", !", geneName))
  }
  
  # Combine clauses with OR
  ruleStr <- paste(clauseList, collapse = " | ")
  return(paste0(geneName, ", ", ruleStr))
}

#' Sanitize gene names inside a Boolean rule string
#'
#' Replaces gene names in a BoolNet-compatible rule string with sanitized
#' versions based on a geneMap data frame. Preserves Boolean logic operators.
#'
#' @param ruleStr Character string in BoolNet format
#' @param geneMap Data frame with columns `originalName` and `sanitizedName`
#' @return Sanitized BoolNet rule string
#' @export
sanitizeRule <- function(ruleStr, geneMap) {
  parts <- strsplit(ruleStr, ",\\s*")[[1]]
  target <- parts[1]
  logic <- parts[2]
  
  # Sort by name length (longest first) to avoid partial matches
  geneNames <- geneMap$OriginalName[order(-nchar(geneMap$OriginalName))]
  
  # Replace gene names in logic, preserving spaces and operators
  for (i in seq_along(geneNames)) {
    origName <- geneNames[i]
    safeName <- geneMap$SanitizedName[geneMap$OriginalName == origName]
    
    # Use word boundary regex to avoid partial matches
    pattern <- paste0("\\b", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", origName), "\\b")
    logic <- gsub(pattern, safeName, logic)
    
    # Also check target
    if (target == origName) {
      target <- safeName
    }
  }
  
  return(paste0(target, ", ", logic))
}

#' Lookup sanitized gene name from a mapping data frame
#'
#' Given an original gene name, returns its sanitized version using a geneMap
#' data frame. If the name is not found, returns the input as-is.
#'
#' @param name Character string (original gene name to sanitize)
#' @param geneMap Data frame with columns `originalName` and `sanitizedName`
#' @return Sanitized gene name (character)
#' @export
lookupSanitized <- function(name, geneMap) {
  matchIdx <- match(name, geneMap$OriginalName)
  if (!is.na(matchIdx)) {
    geneMap$SanitizedName[matchIdx]
  } else {
    name
  }
}

# =============================================================================
# Gene Name Sanitization Utilities
# =============================================================================

#' Sanitize a single gene name for BoolNet compatibility
#'
#' Removes special characters and ensures the name starts with a letter.
#' Used to create BoolNet-compatible gene identifiers.
#'
#' @param name Character string to sanitize
#' @return Sanitized gene name
#' @export
sanitizeGeneName <- function(name) {
  name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  if (grepl("^[0-9]", name2)) name2 <- paste0("X", name2)
  name2
}

#' Generate mapping between original and sanitized gene names
#'
#' Creates a data frame mapping original gene names to sanitized versions
#' suitable for BoolNet networks.
#'
#' @param geneNames Character vector of original gene names
#' @return Data frame with OriginalName and SanitizedName columns
#' @export
generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(
    OriginalName = geneNames, 
    SanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE), 
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Binary State Decoding Utilities
# =============================================================================

#' Decode big integer state to binary vector
#'
#' Converts an encoded state (as big integer) back to a binary vector
#' representing gene ON/OFF states. Used for attractor analysis.
#'
#' @param encodedVal Big integer encoded state
#' @param nGenes Number of genes in the network
#' @return Binary vector of length nGenes
#' @export
decodeBigIntegerState <- function(encodedVal, nGenes) {
  # Note: Requires gmp package for bigz operations
  if (!requireNamespace("gmp", quietly = TRUE)) {
    stop("gmp package required for big integer operations")
  }
  
  valBig <- gmp::as.bigz(encodedVal)
  bitVec <- numeric(nGenes)
  
  for (i in seq_len(nGenes)) {
    bitVec[i] <- as.integer(gmp::mod.bigz(valBig, gmp::as.bigz(2)))
    valBig <- valBig %/% gmp::as.bigz(2)
    if (valBig == 0) break
  }
  
  return(bitVec)
}