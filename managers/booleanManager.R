# =============================================================================
# booleanManager.R
# Purpose: Boolean network construction, rule inference, and optimization
# 
# =============================================================================

# =============================================================================
# Enhanced Boolean Rule Synthesis with Multiple Templates
# =============================================================================

#' Generate multiple Boolean rule templates for a target gene
#' Called by: synthesizeBestBooleanRule
#' KEEP!
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector (1 for activation, -1 for repression)
#' @return List of rule templates with descriptions
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
    andRepressors <- if(length(repressors) > 0) paste0("!", paste(repressors, collapse = " & !")) else ""
    
    templates$balanced_regulation <- list(
      rule = paste0(targetGene, ", ", orActivators, " & ", andRepressors),
      description = "Balanced: any activator + all repressors off"
    )
  }
  
  if (length(activators) >= 2 && length(repressors) >= 1) {
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
#' Called by: synthesizeBooleanRulesBatch
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector of regulatory signs
#' @param matBin Binary expression matrix
#' @return List with best rule, score, and all tested templates
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
#' Called by: Main script(s)!
#' KEEP!
#'
#' @param edges Data frame with columns TF, Target, regType
#' @param matBin Binary expression matrix
#' @param maxRegulators Maximum number of regulators per target
#' @return Named list of rule objects with template information
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
        biologicallyPlausible = TRUE
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
#' Called by: synthesizeBestBooleanRule
#' KEEP!
#'
#' @param ruleStr BoolNet-style rule string
#' @param matBin Binary expression matrix
#' @return Accuracy score (0 to 1)
scoreBooleanRule <- function(ruleStr, matBin) {
  
  # Main scoring logic
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

#' Generate BoolNet-compatible rule string from Boolean function
#'
#' Converts a Boolean function output pattern into a human-readable rule string
#' suitable for BoolNet network specification.
#'
#' @param geneName Name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param outPattern Integer vector representing Boolean function output
#' @return Character string in BoolNet rule format
#' @export
makeBoolNetRule <- function(geneName, regulators, outPattern) {
  
  if (length(regulators) == 0) {
    ruleStr <- ifelse(outPattern[1] == 1, geneName, paste0("! ", geneName))
    return(paste0(geneName, ", ", ruleStr))
  }
  
  R <- length(regulators)
  expectedLength <- 2^R
  
  if (length(outPattern) != expectedLength) {
    stop("Output pattern length (", length(outPattern),
         ") doesn't match expected length (", expectedLength, ") for ", R, " regulators")
  }
  
  clauseList <- character()
  
  for (i in 0:(expectedLength - 1)) {
    if (outPattern[i + 1] == 1) {
      bits <- as.integer(intToBits(i))[1:R]
      andTerms <- ifelse(bits == 1, regulators, paste0("! ", regulators))
      clause <- paste0("(", paste(andTerms, collapse = " & "), ")")
      clauseList <- c(clauseList, clause)
    }
  }
  
  if (length(clauseList) == 0) {
    return(paste0(geneName, ", ! ", geneName))
  }
  
  ruleStr <- paste(clauseList, collapse = " | ")
  return(paste0(geneName, ", ", ruleStr))
}

#' Sanitize gene names inside a Boolean rule string (FIXED VERSION)
#'
#' Replaces gene names in a BoolNet-compatible rule string with sanitized
#' versions based on a geneMap data frame. Preserves Boolean logic operators
#' and proper spacing.
#'
#' @param ruleStr A character string in BoolNet format (e.g., "TP53, !MDM2 & ATM")
#' @param geneMap A data frame with columns `originalName` and `sanitizedName`
#'
#' @return A sanitized BoolNet rule string
sanitizeRule <- function(ruleStr, geneMap) {
  parts <- strsplit(ruleStr, ",\\s*")[[1]]
  target <- parts[1]
  logic  <- parts[2]
  
  # Create a simple lookup for gene name replacement
  # Sort by name length (longest first) to avoid partial matches
  geneNames <- geneMap$originalName[order(-nchar(geneMap$originalName))]
  
  # Replace gene names in logic, preserving spaces and operators
  for (i in seq_along(geneNames)) {
    origName <- geneNames[i]
    safeName <- geneMap$sanitizedName[geneMap$originalName == origName]
    
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
#' data frame with columns `originalName` and `sanitizedName`. If the name is
#' not found in the map, returns the input as-is.
#'
#' @param name A character string (original gene name to sanitize)
#' @param geneMap A data frame with columns `originalName` and `sanitizedName`
#'
#' @return Sanitized gene name (character)
#' @examples
#' lookupSanitized("GATA-1", geneMap)
lookupSanitized <- function(name, geneMap) {
  matchIdx <- match(name, geneMap$originalName)
  if (!is.na(matchIdx)) {
    geneMap$sanitizedName[matchIdx]
  } else {
    name
  }
}
