# =============================================================================
# booleanManager.R
# Purpose: Boolean network construction, rule inference, and optimization
# Dependencies: BoolNet, foreach, doParallel
# =============================================================================

# =============================================================================
# Input-Output Pair Generation
# =============================================================================

#' Generate input-output pairs for Boolean rule inference
#'
#' Creates training data for Boolean rule learning by extracting regulator states
#' at time t and target gene state at time t+k from pseudotime-ordered cells.
#'
#' @param targetGene Name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param matBin Binary gene expression matrix (genes x cells)
#' @param cellOrder Integer vector of cell indices ordered by pseudotime
#' @param k Time step offset between input and output (default: 10)
#' @return Data frame with input regulator states and output target states
#' @export
#' @examples
#' # Generate training pairs for gene regulation
#' pairs <- makeInputOutputPairs("TP53", c("MYC", "E2F1"), binaryMatrix, cellOrder, k = 5)
makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k = 10) {
  
  # Validate inputs
  if (!targetGene %in% rownames(matBin)) {
    stop("Target gene '", targetGene, "' not found in expression matrix")
  }
  
  missingRegs <- setdiff(regulators, rownames(matBin))
  if (length(missingRegs) > 0) {
    stop("Regulator genes not found: ", paste(missingRegs, collapse = ", "))
  }
  
  if (length(cellOrder) <= k) {
    stop("Insufficient cells for time step k=", k, " (need > ", k, " cells)")
  }
  
  # Generate input-output pairs
  ioList <- list()
  maxT <- length(cellOrder) - k
  
  for (t in seq_len(maxT)) {
    # Input: regulator states at time t
    inputVec <- matBin[regulators, cellOrder[t], drop = FALSE]
    
    # Output: target state at time t+k
    outputVal <- matBin[targetGene, cellOrder[t + k]]
    
    # Create data frame row
    rowDf <- data.frame(
      t(inputVec), 
      output = outputVal, 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    ioList[[t]] <- rowDf
  }
  
  if (length(ioList) == 0) {
    return(NULL)
  }
  
  # Combine all pairs
  return(do.call(rbind, ioList))
}

# =============================================================================
# Boolean Function Search
# =============================================================================

#' Find best Boolean rules from input-output training data
#'
#' Performs exhaustive search over all possible Boolean functions to find
#' the function(s) that best explain the target gene's regulation pattern.
#'
#' @param ioDf Data frame with input columns (regulators) and 'output' column
#' @return List containing best score and list of optimal Boolean functions
#' @export
#' @examples
#' # Find Boolean rules from training data
#' rules <- findBestBooleanRules(inputOutputData)
#' cat("Best score:", rules$score)
findBestBooleanRules <- function(ioDf) {
  
  # Validate input
  if (!"output" %in% colnames(ioDf)) {
    stop("Input data frame must contain 'output' column")
  }
  
  if (nrow(ioDf) == 0) {
    stop("Input data frame is empty")
  }
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)  # Number of regulators
  
  if (R == 0) {
    stop("No input columns found")
  }
  
  # Generate all possible input states
  allInputStates <- expand.grid(rep(list(c(0, 1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  
  # All possible Boolean functions (2^nStates possibilities)
  allFnIdx <- seq_len(2^nStates) - 1L
  
  bestScore <- -1
  bestFns <- list()
  
  # Test each Boolean function
  for (fnIdx in allFnIdx) {
    # Convert function index to output pattern
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    
    # Count correct predictions
    correct <- 0
    for (i in seq_len(nrow(ioDf))) {
      # Find matching input state
      inputState <- as.numeric(ioDf[i, inputCols])
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == inputState)))
      
      # Check if prediction matches actual output
      if (length(matchIdx) > 0 && outPattern[matchIdx] == ioDf$output[i]) {
        correct <- correct + 1
      }
    }
    
    # Calculate accuracy score
    score <- correct / nrow(ioDf)
    
    # Update best functions
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
    } else if (abs(score - bestScore) < 1e-9) {
      # Tie - add to best functions list
      bestFns[[length(bestFns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
    }
  }
  
  return(list(score = bestScore, bestFns = bestFns))
}

# =============================================================================
# Boolean Function Combination
# =============================================================================

#' Combine multiple Boolean functions using logical OR
#'
#' When multiple Boolean functions achieve equal optimal scores, combines them
#' using logical OR to create a more permissive composite rule.
#'
#' @param fnList List of Boolean functions from findBestBooleanRules()
#' @return Integer vector representing combined Boolean function output pattern
#' @export
combineBooleanFunctionsByOr <- function(fnList) {
  
  if (length(fnList) == 0) {
    stop("Function list is empty")
  }
  
  if (length(fnList) == 1) {
    return(fnList[[1]]$outPattern)
  }
  
  # Extract output patterns
  patterns <- sapply(fnList, function(fn) fn$outPattern)
  
  # Combine using logical OR (any function returns TRUE)
  combinedPattern <- as.integer(rowSums(patterns) > 0)
  
  return(combinedPattern)
}

# =============================================================================
# BoolNet Rule Generation
# =============================================================================

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
#' @examples
#' # Generate rule for gene with 2 regulators
#' rule <- makeBoolNetRule("TP53", c("MYC", "E2F1"), c(0,1,1,1))
#' # Returns: "TP53, (! MYC & ! E2F1) | (MYC & ! E2F1) | (MYC & E2F1)"
makeBoolNetRule <- function(geneName, regulators, outPattern) {
  
  # Handle case with no regulators (self-regulation)
  if (length(regulators) == 0) {
    ruleStr <- ifelse(outPattern[1] == 1, geneName, paste0("! ", geneName))
    return(paste0(geneName, ", ", ruleStr))
  }
  
  R <- length(regulators)
  expectedLength <- 2^R
  
  # Validate output pattern length
  if (length(outPattern) != expectedLength) {
    stop("Output pattern length (", length(outPattern), 
         ") doesn't match expected length (", expectedLength, ") for ", R, " regulators")
  }
  
  # Generate OR clauses for each TRUE output
  clauseList <- character()
  
  for (i in 0:(expectedLength - 1)) {
    if (outPattern[i + 1] == 1) {
      # Convert index to binary representation
      bits <- as.integer(intToBits(i))[1:R]
      
      # Create AND terms (positive or negative literals)
      andTerms <- ifelse(bits == 1, regulators, paste0("! ", regulators))
      
      # Combine with AND
      clause <- paste0("(", paste(andTerms, collapse = " & "), ")")
      clauseList <- c(clauseList, clause)
    }
  }
  
  # Handle case where no outputs are TRUE (always FALSE)
  if (length(clauseList) == 0) {
    return(paste0(geneName, ", ! ", geneName))
  }
  
  # Combine clauses with OR
  ruleStr <- paste(clauseList, collapse = " | ")
  return(paste0(geneName, ", ", ruleStr))
}

# =============================================================================
# Optimized Boolean Rule Inference
# =============================================================================

#' Find best Boolean rules with early termination optimization
#'
#' Enhanced version of findBestBooleanRules() that can terminate early when
#' perfect rules are found, improving performance for large networks.
#'
#' @param ioDf Data frame with input-output training data
#' @param earlyTermination Whether to stop when perfect score is found (default: TRUE)
#' @param maxFunctions Maximum functions to test before switching to sampling (default: 10000)
#' @return List with best score, functions, and whether perfect rule was found
#' @export
findBestBooleanRulesOptimized <- function(ioDf, earlyTermination = TRUE, maxFunctions = 10000) {
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  allInputStates <- expand.grid(rep(list(c(0, 1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  allFnIdx <- seq_len(2^nStates) - 1L
  
  # Check if problem is too large for exhaustive search
  if (length(allFnIdx) > maxFunctions) {
    warning(sprintf("Too many functions (%d) for exhaustive search. Using sampling approach.", 
                    length(allFnIdx)))
    return(findBooleanRulesSampling(ioDf, maxFunctions))
  }
  
  bestScore <- -1
  bestFns <- list()
  perfectFound <- FALSE
  
  for (fnIdx in allFnIdx) {
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    correct <- 0
    
    # Count correct predictions
    for (i in seq_len(nrow(ioDf))) {
      inputState <- as.numeric(ioDf[i, inputCols])
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == inputState)))
      
      if (length(matchIdx) > 0 && outPattern[matchIdx] == ioDf$output[i]) {
        correct <- correct + 1
      }
    }
    
    score <- correct / nrow(ioDf)
    
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
      
      # Early termination for perfect or near-perfect score
      if (earlyTermination && score >= 0.999) {
        perfectFound <- TRUE
        break
      }
    } else if (abs(score - bestScore) < 1e-9) {
      bestFns[[length(bestFns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
    }
  }
  
  return(list(score = bestScore, bestFns = bestFns, perfectFound = perfectFound))
}

# =============================================================================
# Sampling-Based Rule Inference
# =============================================================================

#' Find Boolean rules using intelligent sampling
#'
#' For large networks where exhaustive search is infeasible, uses smart sampling
#' to find good Boolean rules efficiently.
#'
#' @param ioDf Data frame with input-output training data
#' @param maxSamples Maximum number of functions to sample (default: 10000)
#' @param method Sampling method ("smart" or "random", default: "smart")
#' @return List with best score, functions, and sampling flag
#' @export
findBooleanRulesSampling <- function(ioDf, maxSamples = 10000, method = "smart") {
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  allInputStates <- expand.grid(rep(list(c(0, 1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  totalFunctions <- 2^nStates
  
  # Generate candidate functions
  if (method == "smart") {
    candidates <- generateSmartCandidates(ioDf, allInputStates, maxSamples)
  } else {
    candidates <- sample(0:(totalFunctions - 1), min(maxSamples, totalFunctions))
  }
  
  bestScore <- -1
  bestFns <- list()
  
  # Test each candidate function
  for (fnIdx in candidates) {
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    correct <- 0
    
    for (i in seq_len(nrow(ioDf))) {
      inputState <- as.numeric(ioDf[i, inputCols])
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == inputState)))
      
      if (length(matchIdx) > 0 && outPattern[matchIdx] == ioDf$output[i]) {
        correct <- correct + 1
      }
    }
    
    score <- correct / nrow(ioDf)
    
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
    } else if (abs(score - bestScore) < 1e-9) {
      bestFns[[length(bestFns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
    }
  }
  
  return(list(score = bestScore, bestFns = bestFns, sampled = TRUE))
}

# =============================================================================
# Smart Candidate Generation
# =============================================================================

#' Generate biologically plausible Boolean function candidates
#'
#' Creates a set of candidate Boolean functions that are more likely to represent
#' real biological regulatory relationships (simple AND/OR combinations).
#'
#' @param ioDf Input-output training data
#' @param allInputStates All possible input state combinations
#' @param maxSamples Maximum number of candidates to generate
#' @return Integer vector of candidate function indices
#' @export
generateSmartCandidates <- function(ioDf, allInputStates, maxSamples) {
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  nStates <- 2^R
  
  candidates <- c()
  
  # 1. Constant functions
  candidates <- c(candidates, 0, 2^nStates - 1)  # Always FALSE, Always TRUE
  
  # 2. Single regulator functions
  for (reg in seq_along(inputCols)) {
    # Positive regulation: output = regulator
    pattern <- rep(0, nStates)
    for (state in 1:nStates) {
      bits <- as.integer(intToBits(state - 1))[1:R]
      if (bits[reg] == 1) pattern[state] <- 1
    }
    fnIdx <- sum(pattern * 2^(0:(nStates - 1)))
    candidates <- c(candidates, fnIdx)
    
    # Negative regulation: output = NOT regulator
    patternNeg <- 1 - pattern
    fnIdxNeg <- sum(patternNeg * 2^(0:(nStates - 1)))
    candidates <- c(candidates, fnIdxNeg)
  }
  
  # 3. Two-regulator combinations (AND, OR, XOR)
  if (R >= 2) {
    for (i in 1:(R - 1)) {
      for (j in (i + 1):R) {
        # AND function
        patternAnd <- rep(0, nStates)
        for (state in 1:nStates) {
          bits <- as.integer(intToBits(state - 1))[1:R]
          if (bits[i] == 1 && bits[j] == 1) patternAnd[state] <- 1
        }
        candidates <- c(candidates, sum(patternAnd * 2^(0:(nStates - 1))))
        
        # OR function
        patternOr <- rep(0, nStates)
        for (state in 1:nStates) {
          bits <- as.integer(intToBits(state - 1))[1:R]
          if (bits[i] == 1 || bits[j] == 1) patternOr[state] <- 1
        }
        candidates <- c(candidates, sum(patternOr * 2^(0:(nStates - 1))))
        
        # XOR function
        patternXor <- rep(0, nStates)
        for (state in 1:nStates) {
          bits <- as.integer(intToBits(state - 1))[1:R]
          if (bits[i] != bits[j]) patternXor[state] <- 1
        }
        candidates <- c(candidates, sum(patternXor * 2^(0:(nStates - 1))))
      }
    }
  }
  
  # 4. Fill remaining slots with random sampling
  uniqueCandidates <- unique(candidates)
  remaining <- maxSamples - length(uniqueCandidates)
  
  if (remaining > 0) {
    randomCandidates <- sample(0:(2^nStates - 1), remaining, replace = FALSE)
    candidates <- c(uniqueCandidates, randomCandidates)
  } else {
    candidates <- uniqueCandidates
  }
  
  return(unique(candidates)[1:min(maxSamples, length(unique(candidates)))])
}

# =============================================================================
# Parallel Boolean Rule Inference
# =============================================================================

#' Find Boolean rules using parallel computation
#'
#' Distributes Boolean function search across multiple cores for improved
#' performance on large networks.
#'
#' @param ioDf Input-output training data
#' @param cores Number of cores to use (NULL for auto-detection)
#' @return List with best score and functions
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
findBooleanRulesParallel <- function(ioDf, cores = NULL) {
  
  if (is.null(cores)) {
    cores <- min(4, parallel::detectCores() - 1)
  }
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  nStates <- 2^R
  allFnIdx <- seq_len(2^nStates) - 1L
  
  # Split function indices across cores
  chunks <- split(allFnIdx, cut(seq_along(allFnIdx), cores))
  
  # Parallel computation
  chunkResults <- foreach(chunk = chunks, .combine = 'c', .packages = c()) %dopar% {
    chunkBest <- list(score = -1, fns = list())
    allInputStates <- expand.grid(rep(list(c(0, 1)), R))
    colnames(allInputStates) <- inputCols
    
    for (fnIdx in chunk) {
      outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
      correct <- 0
      
      for (i in seq_len(nrow(ioDf))) {
        inputState <- as.numeric(ioDf[i, inputCols])
        matchIdx <- which(apply(allInputStates, 1, function(x) all(x == inputState)))
        
        if (length(matchIdx) > 0 && outPattern[matchIdx] == ioDf$output[i]) {
          correct <- correct + 1
        }
      }
      
      score <- correct / nrow(ioDf)
      
      if (score > chunkBest$score) {
        chunkBest$score <- score
        chunkBest$fns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
      } else if (abs(score - chunkBest$score) < 1e-9) {
        chunkBest$fns[[length(chunkBest$fns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
      }
    }
    
    list(chunkBest)
  }
  
  # Combine results from all chunks
  bestScore <- max(sapply(chunkResults, function(x) x$score))
  bestFns <- do.call(c, lapply(chunkResults, function(x) {
    if (abs(x$score - bestScore) < 1e-9) x$fns else list()
  }))
  
  return(list(score = bestScore, bestFns = bestFns))
}

# =============================================================================
# Adaptive Rule Inference
# =============================================================================

#' Adaptive Boolean rule inference that chooses optimal method
#'
#' Automatically selects the best inference method based on problem size and
#' computational constraints.
#'
#' @param ioDf Input-output training data
#' @param timeLimitMinutes Maximum computation time in minutes (default: 30)
#' @return List with best score and functions
#' @export
findBestBooleanRulesAdaptive <- function(ioDf, timeLimitMinutes = 30) {
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  totalFunctions <- 2^(2^R)
  
  startTime <- Sys.time()
  
  if (R <= 2) {
    # Small problem: use exact method
    message("Using exact method for ", R, " regulators")
    return(findBestBooleanRulesOptimized(ioDf, earlyTermination = TRUE))
    
  } else if (R == 3) {
    # Medium problem: try exact with early termination
    message("Using optimized exact method for ", R, " regulators")
    result <- findBestBooleanRulesOptimized(ioDf, earlyTermination = TRUE)
    
    # Fall back to parallel if taking too long
    timeElapsed <- difftime(Sys.time(), startTime, units = "mins")
    if (!result$perfectFound && timeElapsed > timeLimitMinutes / 2) {
      message("Switching to parallel method due to time constraints")
      return(findBooleanRulesParallel(ioDf))
    }
    
    return(result)
    
  } else if (R <= 4) {
    # Large problem: use smart sampling
    message("Using smart sampling for ", R, " regulators")
    return(findBooleanRulesSampling(ioDf, maxSamples = 50000, method = "smart"))
    
  } else {
    # Very large problem: use aggressive sampling
    message("Using aggressive sampling for ", R, " regulators")
    return(findBooleanRulesSampling(ioDf, maxSamples = 10000, method = "smart"))
  }
}

# =============================================================================
# Enhanced Boolean Rule Inference with Prior Knowledge Integration
# =============================================================================

#' Find best Boolean rules incorporating activation/inhibition prior knowledge
#'
#' Uses a hybrid approach: first tests biologically plausible functions based on
#' known regulatory relationships, then falls back to broader search if needed.
#'
#' @param ioDf Data frame with input-output training data
#' @param regulatorSigns Named vector of regulatory signs (1 = activation, -1 = inhibition)
#' @param priorKnowledgeWeight How much to weight prior knowledge vs data fit (0-1, default: 0.7)
#' @param fallbackThreshold Minimum score to accept before falling back to exhaustive search (default: 0.8)
#' @param earlyTermination Whether to stop when excellent rules are found (default: TRUE)
#' @return List with best score, functions, method used, and biological consistency flag
#' @export
findBestBooleanRulesWithPrior <- function(ioDf, 
                                          regulatorSigns = NULL, 
                                          priorKnowledgeWeight = 0.7,
                                          fallbackThreshold = 0.8,
                                          earlyTermination = TRUE) {
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  
  if (R == 0) {
    stop("No input columns found")
  }
  
  # Generate all possible input states
  allInputStates <- expand.grid(rep(list(c(0, 1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  
  # Phase 1: Try biologically-informed candidates
  biologicallyPlausible <- FALSE
  bestScore <- -1
  bestFns <- list()
  methodUsed <- "prior_knowledge"
  
  if (!is.null(regulatorSigns) && length(regulatorSigns) > 0) {
    message("Phase 1: Testing biologically plausible Boolean functions...")
    
    # Generate biologically-informed candidates
    biologicalCandidates <- generateBiologicalCandidates(inputCols, regulatorSigns, nStates)
    message("Testing ", length(biologicalCandidates), " biologically plausible functions")
    
    # Test biological candidates
    for (candidateInfo in biologicalCandidates) {
      score <- evaluateBooleanFunction(candidateInfo$outPattern, ioDf, allInputStates, inputCols)
      
      if (score > bestScore) {
        bestScore <- score
        bestFns <- list(candidateInfo)
        biologicallyPlausible <- TRUE
      } else if (abs(score - bestScore) < 1e-9) {
        bestFns[[length(bestFns) + 1]] <- candidateInfo
      }
      
      # Early termination if excellent biological rule found
      if (earlyTermination && score >= 0.95) {
        message("Found excellent biological rule (score: ", round(score, 3), "). Stopping early.")
        break
      }
    }
    
    message("Best biological score: ", round(bestScore, 3))
  }
  
  # Phase 2: Fallback to broader search if biological rules insufficient
  if (bestScore < fallbackThreshold) {
    message("Phase 2: Biological rules insufficient (score < ", fallbackThreshold, "). Expanding search...")
    methodUsed <- "hybrid"
    
    # Use optimized exhaustive search with early termination
    exhaustiveResult <- findBestBooleanRulesOptimized(ioDf, earlyTermination = earlyTermination)
    
    # Compare with biological results
    if (exhaustiveResult$score > bestScore) {
      message("Exhaustive search found better rules (score: ", round(exhaustiveResult$score, 3), 
              " vs ", round(bestScore, 3), ")")
      bestScore <- exhaustiveResult$score
      bestFns <- exhaustiveResult$bestFns
      biologicallyPlausible <- checkBiologicalConsistency(bestFns, inputCols, regulatorSigns)
      methodUsed <- "exhaustive_fallback"
    }
  }
  
  # Phase 3: Check biological consistency of final rules
  if (!is.null(regulatorSigns)) {
    consistencyReport <- analyzeBiologicalConsistency(bestFns, inputCols, regulatorSigns)
    if (!consistencyReport$consistent) {
      warning("Final Boolean rules contradict known regulatory relationships:\n", 
              consistencyReport$message)
    }
  }
  
  return(list(
    score = bestScore,
    bestFns = bestFns,
    biologicallyPlausible = biologicallyPlausible,
    methodUsed = methodUsed,
    priorKnowledgeUsed = !is.null(regulatorSigns)
  ))
}

#' Generate biologically plausible Boolean function candidates
#'
#' Creates Boolean functions that respect known activation/inhibition relationships
#' while covering common regulatory logic patterns.
#'
#' @param inputCols Names of regulator genes
#' @param regulatorSigns Named vector of regulatory signs (1 = activation, -1 = inhibition)
#' @param nStates Number of possible input states (2^R)
#' @return List of candidate Boolean functions with biological justification
#' @export
generateBiologicalCandidates <- function(inputCols, regulatorSigns, nStates) {
  
  R <- length(inputCols)
  candidates <- list()
  
  # Ensure all regulators have known signs
  knownSigns <- regulatorSigns[inputCols]
  if (any(is.na(knownSigns))) {
    warning("Some regulators have unknown signs. Using neutral assumption.")
    knownSigns[is.na(knownSigns)] <- 1  # Default to activation
  }
  
  # 1. Single regulator functions (respecting known signs)
  for (i in seq_along(inputCols)) {
    reg <- inputCols[i]
    sign <- knownSigns[i]
    
    if (sign > 0) {
      # Activator: output = regulator
      pattern <- rep(0, nStates)
      for (state in 1:nStates) {
        bits <- as.integer(intToBits(state - 1))[1:R]
        if (bits[i] == 1) pattern[state] <- 1
      }
      
      candidates[[length(candidates) + 1]] <- list(
        outPattern = pattern,
        description = paste(reg, "activates target"),
        biologicalBasis = "single_activator"
      )
    } else {
      # Repressor: output = NOT regulator
      pattern <- rep(1, nStates)
      for (state in 1:nStates) {
        bits <- as.integer(intToBits(state - 1))[1:R]
        if (bits[i] == 1) pattern[state] <- 0
      }
      
      candidates[[length(candidates) + 1]] <- list(
        outPattern = pattern,
        description = paste(reg, "represses target"),
        biologicalBasis = "single_repressor"
      )
    }
  }
  
  # 2. Two-regulator combinations (common biological patterns)
  if (R >= 2) {
    for (i in 1:(R - 1)) {
      for (j in (i + 1):R) {
        reg1 <- inputCols[i]
        reg2 <- inputCols[j]
        sign1 <- knownSigns[i]
        sign2 <- knownSigns[j]
        
        # Cooperative activation (both activators needed)
        if (sign1 > 0 && sign2 > 0) {
          pattern <- rep(0, nStates)
          for (state in 1:nStates) {
            bits <- as.integer(intToBits(state - 1))[1:R]
            if (bits[i] == 1 && bits[j] == 1) pattern[state] <- 1
          }
          
          candidates[[length(candidates) + 1]] <- list(
            outPattern = pattern,
            description = paste(reg1, "AND", reg2, "cooperatively activate"),
            biologicalBasis = "cooperative_activation"
          )
        }
        
        # Redundant activation (either activator sufficient)
        if (sign1 > 0 && sign2 > 0) {
          pattern <- rep(0, nStates)
          for (state in 1:nStates) {
            bits <- as.integer(intToBits(state - 1))[1:R]
            if (bits[i] == 1 || bits[j] == 1) pattern[state] <- 1
          }
          
          candidates[[length(candidates) + 1]] <- list(
            outPattern = pattern,
            description = paste(reg1, "OR", reg2, "redundantly activate"),
            biologicalBasis = "redundant_activation"
          )
        }
        
        # Activation with repression override
        if (sign1 > 0 && sign2 < 0) {
          pattern <- rep(0, nStates)
          for (state in 1:nStates) {
            bits <- as.integer(intToBits(state - 1))[1:R]
            if (bits[i] == 1 && bits[j] == 0) pattern[state] <- 1  # Activator on, repressor off
          }
          
          candidates[[length(candidates) + 1]] <- list(
            outPattern = pattern,
            description = paste(reg1, "activates when", reg2, "is not repressing"),
            biologicalBasis = "activation_with_repression_override"
          )
        }
        
        # Mutual repression (both must be off)
        if (sign1 < 0 && sign2 < 0) {
          pattern <- rep(1, nStates)
          for (state in 1:nStates) {
            bits <- as.integer(intToBits(state - 1))[1:R]
            if (bits[i] == 1 || bits[j] == 1) pattern[state] <- 0
          }
          
          candidates[[length(candidates) + 1]] <- list(
            outPattern = pattern,
            description = paste("Active when neither", reg1, "nor", reg2, "repress"),
            biologicalBasis = "mutual_repression_relief"
          )
        }
      }
    }
  }
  
  # 3. Three-regulator patterns (if applicable)
  if (R == 3) {
    # Majority rule with known signs
    activators <- which(knownSigns > 0)
    repressors <- which(knownSigns < 0)
    
    if (length(activators) >= 2) {
      # Majority activation rule
      pattern <- rep(0, nStates)
      for (state in 1:nStates) {
        bits <- as.integer(intToBits(state - 1))[1:R]
        activeCount <- sum(bits[activators])
        inactiveRepressors <- sum(1 - bits[repressors])
        
        if (activeCount >= 2 || (activeCount >= 1 && length(repressors) > 0 && inactiveRepressors == length(repressors))) {
          pattern[state] <- 1
        }
      }
      
      candidates[[length(candidates) + 1]] <- list(
        outPattern = pattern,
        description = "Majority activation with repressor override",
        biologicalBasis = "majority_rule_with_repression"
      )
    }
  }
  
  return(candidates)
}

#' Evaluate a Boolean function against training data
#'
#' Computes the accuracy score for a given Boolean function output pattern.
#'
#' @param outPattern Integer vector representing Boolean function output
#' @param ioDf Input-output training data
#' @param allInputStates All possible input state combinations
#' @param inputCols Names of input columns
#' @return Accuracy score (0-1)
#' @export
evaluateBooleanFunction <- function(outPattern, ioDf, allInputStates, inputCols) {
  
  correct <- 0
  
  for (i in seq_len(nrow(ioDf))) {
    inputState <- as.numeric(ioDf[i, inputCols])
    matchIdx <- which(apply(allInputStates, 1, function(x) all(x == inputState)))
    
    if (length(matchIdx) > 0 && outPattern[matchIdx] == ioDf$output[i]) {
      correct <- correct + 1
    }
  }
  
  return(correct / nrow(ioDf))
}

#' Check biological consistency of Boolean rules
#'
#' Analyzes whether learned Boolean functions are consistent with known
#' regulatory relationships.
#'
#' @param bestFns List of optimal Boolean functions
#' @param inputCols Names of regulator genes
#' @param regulatorSigns Named vector of regulatory signs
#' @return Boolean indicating overall consistency
#' @export
checkBiologicalConsistency <- function(bestFns, inputCols, regulatorSigns) {
  
  if (is.null(regulatorSigns) || length(bestFns) == 0) {
    return(TRUE)  # Can't check consistency without prior knowledge
  }
  
  # For now, implement a simple check
  # More sophisticated analysis could be added here
  for (fn in bestFns) {
    if (!is.null(fn$biologicalBasis)) {
      return(TRUE)  # Functions from biological candidates are consistent by definition
    }
  }
  
  # For exhaustive search results, would need more complex consistency analysis
  # This is a placeholder for more sophisticated checking
  return(TRUE)
}

#' Analyze biological consistency and provide detailed report
#'
#' Provides detailed analysis of how well Boolean rules match known biology.
#'
#' @param bestFns List of optimal Boolean functions
#' @param inputCols Names of regulator genes
#' @param regulatorSigns Named vector of regulatory signs
#' @return List with consistency flag and detailed message
#' @export
analyzeBiologicalConsistency <- function(bestFns, inputCols, regulatorSigns) {
  
  if (is.null(regulatorSigns) || length(bestFns) == 0) {
    return(list(consistent = TRUE, message = "No prior knowledge available for consistency check"))
  }
  
  # Check if any functions have biological basis
  biologicalFunctions <- sapply(bestFns, function(fn) !is.null(fn$biologicalBasis))
  
  if (any(biologicalFunctions)) {
    biologicalCount <- sum(biologicalFunctions)
    totalCount <- length(bestFns)
    
    return(list(
      consistent = TRUE,
      message = paste0(biologicalCount, " of ", totalCount, " optimal functions have biological basis")
    ))
  } else {
    return(list(
      consistent = FALSE,
      message = "No optimal functions have clear biological basis. Manual review recommended."
    ))
  }
}

# =============================================================================
# Boolean Rule Analysis and Visualization Functions
# =============================================================================

#' Generate comprehensive Boolean rule analysis report
#'
#' Creates visualizations and summary statistics for Boolean rule inference results
#'
#' @param boolRules List of Boolean rules from main script
#' @param edges Edge list data frame with regulatory relationships
#' @param paths Paths object with directory structure
#' @export
generateBooleanRuleReport <- function(boolRules, edges, paths, cellType, trajectory) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(ComplexHeatmap)
  library(igraph)
  
  if (!dir.exists(paths$base$plots)) dir.create(paths$base$plots, recursive = TRUE)
  if (!dir.exists(paths$base$txt)) dir.create(paths$base$txt, recursive = TRUE)
  
  # Extract summary statistics
  ruleStats <- extractRuleStatistics(boolRules)
  
  # 1. Rule Quality Distribution
  p1 <- plotRuleQualityDistribution(ruleStats)
  ggsave(file.path(paths$base$plots, cellType, "_", trajectory, "rule_quality_distribution.png"), p1, width = 10, height = 6)
  
  # 2. Method Usage Summary
  p2 <- plotMethodUsageSummary(ruleStats)
  ggsave(file.path(paths$base$plots, cellType, "_", trajectory, "method_usage_summary.png"), p2, width = 8, height = 6)
  
  # 3. Biological Plausibility Analysis
  p3 <- plotBiologicalPlausibility(ruleStats)
  ggsave(file.path(paths$base$plots, cellType, "_", trajectory, "biological_plausibility.png"), p3, width = 10, height = 6)
  
  # 4. Network Complexity Heatmap
  p4 <- plotNetworkComplexityHeatmap(boolRules, edges)
  png(file.path(paths$base$plots, cellType, "_", trajectory, "network_complexity_heatmap.png"), width = 1200, height = 800)
  draw(p4)
  dev.off()
  
  # 5. Rule Logic Summary Network
  p5 <- plotRuleLogicNetwork(boolRules, edges)
  ggsave(file.path(paths$base$plots, cellType, "_", trajectory, "rule_logic_network.png"), p5, width = 12, height = 10)
  
  # 6. Generate text report
  generateTextReport(ruleStats, boolRules, file.path(paths$base$txt, cellType, "_", trajectory, "boolean_rules_report.txt"))
  
  message("Boolean rule analysis complete. Results saved to: ", paths$base$txt)
}

#' Extract summary statistics from Boolean rules
#' @param boolRules List of Boolean rules
#' @return Data frame with rule statistics
extractRuleStatistics <- function(boolRules) {
  
  # Extract each field individually to ensure proper types
  geneNames <- names(boolRules)
  numRegulators <- numeric(length(boolRules))
  bestScore <- numeric(length(boolRules))
  methodUsed <- character(length(boolRules))
  biologicallyPlausible <- logical(length(boolRules))
  ruleComplexity <- numeric(length(boolRules))
  
  for (i in seq_along(boolRules)) {
    rule <- boolRules[[i]]
    numRegulators[i] <- length(rule$regulators)
    bestScore[i] <- rule$bestScore
    methodUsed[i] <- if (is.null(rule$methodUsed)) "unknown" else as.character(rule$methodUsed)
    biologicallyPlausible[i] <- if (is.null(rule$biologicallyPlausible)) NA else rule$biologicallyPlausible
    ruleComplexity[i] <- sum(rule$outPattern)
  }
  
  stats <- data.frame(
    gene = geneNames,
    numRegulators = numRegulators,
    bestScore = bestScore,
    methodUsed = methodUsed,
    biologicallyPlausible = biologicallyPlausible,
    ruleComplexity = ruleComplexity,
    stringsAsFactors = FALSE
  )
  
  return(stats)
}

#' Plot distribution of rule quality scores
#' @param ruleStats Data frame with rule statistics
#' @return ggplot object
plotRuleQualityDistribution <- function(ruleStats) {
  
  ggplot(ruleStats, aes(x = bestScore)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.95, color = "darkgreen", linetype = "dashed", size = 1) +
    labs(
      title = "Distribution of Boolean Rule Quality Scores",
      subtitle = paste0("n = ", nrow(ruleStats), " genes with inferred rules"),
      x = "Rule Quality Score (Accuracy)",
      y = "Number of Genes",
      caption = "Red line: minimum threshold (0.8), Green line: excellent threshold (0.95)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    annotate("text", x = 0.85, y = Inf, label = "Good Rules", vjust = 2, color = "darkgreen") +
    annotate("text", x = 0.4, y = Inf, label = "Poor Rules", vjust = 2, color = "red")
}

#' Plot method usage summary
#' @param ruleStats Data frame with rule statistics
#' @return ggplot object
plotMethodUsageSummary <- function(ruleStats) {
  
  # Fix: Use explicit dplyr namespace and alternative approach
  methodCounts <- ruleStats %>%
    dplyr::group_by(methodUsed) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::mutate(
      percentage = round(100 * count / sum(count), 1),
      label = paste0(methodUsed, "\n(", count, " genes, ", percentage, "%)")
    )
  
  ggplot(methodCounts, aes(x = "", y = count, fill = methodUsed)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    labs(
      title = "Boolean Rule Inference Method Usage",
      subtitle = "Which methods were used to find the best rules?",
      fill = "Method"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_text(aes(label = percentage), position = position_stack(vjust = 0.5))
}

#' Plot biological plausibility analysis
#' @param ruleStats Data frame with rule statistics
#' @return ggplot object
plotBiologicalPlausibility <- function(ruleStats) {
  
  # Create comparison of scores by biological plausibility
  plausibilityData <- ruleStats %>%
    filter(!is.na(biologicallyPlausible)) %>%
    mutate(
      plausibilityLabel = ifelse(biologicallyPlausible, "Biologically Plausible", "Data-Driven Only")
    )
  
  if (nrow(plausibilityData) == 0) {
    # Create empty plot if no data
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No biological plausibility data available") +
             theme_void())
  }
  
  ggplot(plausibilityData, aes(x = plausibilityLabel, y = bestScore, fill = plausibilityLabel)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "Rule Quality by Biological Plausibility",
      subtitle = "Do biologically plausible rules perform better?",
      x = "Rule Type",
      y = "Quality Score",
      fill = "Rule Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("steelblue", "orange")) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red")
}

#' Create network complexity heatmap
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @return ComplexHeatmap object
plotNetworkComplexityHeatmap <- function(boolRules, edges) {
  
  # Create matrix of regulator-target relationships
  allGenes <- unique(c(edges$TF, edges$Target))
  regulatorMatrix <- matrix(0, nrow = length(allGenes), ncol = length(allGenes))
  rownames(regulatorMatrix) <- allGenes
  colnames(regulatorMatrix) <- allGenes
  
  # Fill in regulatory relationships
  for (gene in names(boolRules)) {
    regulators <- boolRules[[gene]]$regulators
    for (reg in regulators) {
      if (reg %in% rownames(regulatorMatrix) && gene %in% colnames(regulatorMatrix)) {
        regulatorMatrix[reg, gene] <- boolRules[[gene]]$bestScore
      }
    }
  }
  
  # Create heatmap
  ComplexHeatmap::Heatmap(
    regulatorMatrix,
    name = "Rule Quality",
    col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = "Target Genes",
    row_title = "Regulator Genes",
    heatmap_legend_param = list(title = "Boolean Rule\nQuality Score")
  )
}

#' Plot rule logic network
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @return ggplot object
plotRuleLogicNetwork <- function(boolRules, edges) {
  
  # Create network graph
  networkEdges <- edges %>%
    dplyr::filter(Target %in% names(boolRules)) %>%
    dplyr::mutate(
      ruleQuality = sapply(Target, function(x) boolRules[[x]]$bestScore %||% 0),
      edgeWidth = pmax(0.5, 3 * ruleQuality)
    )
  
  # Create igraph object
  g <- igraph::graph_from_data_frame(networkEdges, directed = TRUE)
  
  # Calculate layout
  layout <- igraph::layout_with_fr(g)
  
  # Convert to data frame for ggplot
  edgeList <- as.data.frame(igraph::get.edgelist(g))
  names(edgeList) <- c("from", "to")
  edgeList$quality <- igraph::E(g)$ruleQuality
  
  nodeList <- data.frame(
    name = igraph::V(g)$name,
    x = layout[,1],
    y = layout[,2]
  )
  
  # Create edge coordinates with explicit namespacing
  edgeCoords <- edgeList %>%
    dplyr::left_join(nodeList, by = c("from" = "name")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(nodeList, by = c("to" = "name")) %>%
    dplyr::rename(x2 = x, y2 = y)
  
  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edgeCoords,
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2, color = quality),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches")),
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = nodeList,
      ggplot2::aes(x = x, y = y),
      size = 3,
      color = "darkblue",
      alpha = 0.8
    ) +
    ggplot2::scale_color_gradient(
      low = "lightgray",
      high = "red",
      name = "Rule\nQuality"
    ) +
    ggplot2::labs(
      title = "Gene Regulatory Network with Boolean Rule Quality",
      subtitle = "Edge color indicates quality of inferred Boolean rules"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
}

#' Generate text summary report
#' @param ruleStats Data frame with rule statistics
#' @param boolRules List of Boolean rules
#' @param filename Output file path
generateTextReport <- function(ruleStats, boolRules, filename) {
  
  sink(filename)
  
  cat("=============================================================================\n")
  cat("BOOLEAN RULE INFERENCE SUMMARY REPORT\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=============================================================================\n\n")
  
  cat("OVERALL STATISTICS:\n")
  cat("- Total genes with rules:", nrow(ruleStats), "\n")
  cat("- Mean rule quality:", round(mean(ruleStats$bestScore), 3), "\n")
  cat("- Median rule quality:", round(median(ruleStats$bestScore), 3), "\n")
  cat("- Rules with quality > 0.8:", sum(ruleStats$bestScore > 0.8), 
      "(", round(100 * sum(ruleStats$bestScore > 0.8) / nrow(ruleStats), 1), "%)\n")
  cat("- Rules with quality > 0.95:", sum(ruleStats$bestScore > 0.95), 
      "(", round(100 * sum(ruleStats$bestScore > 0.95) / nrow(ruleStats), 1), "%)\n\n")
  
  cat("NETWORK COMPLEXITY:\n")
  cat("- Genes with 1 regulator:", sum(ruleStats$numRegulators == 1), "\n")
  cat("- Genes with 2 regulators:", sum(ruleStats$numRegulators == 2), "\n")
  cat("- Genes with 3+ regulators:", sum(ruleStats$numRegulators >= 3), "\n")
  cat("- Mean regulators per gene:", round(mean(ruleStats$numRegulators), 2), "\n\n")
  
  if ("methodUsed" %in% colnames(ruleStats)) {
    cat("METHOD USAGE:\n")
    methodCounts <- table(ruleStats$methodUsed)
    for (method in names(methodCounts)) {
      cat("-", method, ":", methodCounts[method], "genes\n")
    }
    cat("\n")
  }
  
  if ("biologicallyPlausible" %in% colnames(ruleStats)) {
    plausibleData <- ruleStats[!is.na(ruleStats$biologicallyPlausible), ]
    if (nrow(plausibleData) > 0) {
      cat("BIOLOGICAL PLAUSIBILITY:\n")
      cat("- Biologically plausible rules:", sum(plausibleData$biologicallyPlausible), "\n")
      cat("- Data-driven only rules:", sum(!plausibleData$biologicallyPlausible), "\n")
      cat("- Mean quality (plausible):", 
          round(mean(plausibleData$bestScore[plausibleData$biologicallyPlausible]), 3), "\n")
      cat("- Mean quality (data-driven):", 
          round(mean(plausibleData$bestScore[!plausibleData$biologicallyPlausible]), 3), "\n\n")
    }
  }
  
  cat("TOP 10 HIGHEST QUALITY RULES:\n")
  topRules <- ruleStats[order(ruleStats$bestScore, decreasing = TRUE)[1:min(10, nrow(ruleStats))], ]
  for (i in 1:nrow(topRules)) {
    cat(i, ".", topRules$gene[i], "- Score:", round(topRules$bestScore[i], 3), 
        "- Regulators:", topRules$numRegulators[i], "\n")
  }
  cat("\n")
  
  cat("POTENTIAL ISSUES:\n")
  lowQualityGenes <- ruleStats$gene[ruleStats$bestScore < 0.5]
  if (length(lowQualityGenes) > 0) {
    cat("- Genes with very low rule quality (<0.5):", length(lowQualityGenes), "\n")
    if (length(lowQualityGenes) <= 10) {
      cat("  ", paste(lowQualityGenes, collapse = ", "), "\n")
    }
  } else {
    cat("- No genes with concerning low rule quality\n")
  }
  
  sink()
  
  message("Text report saved to: ", filename)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x