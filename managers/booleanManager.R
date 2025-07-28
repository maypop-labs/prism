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