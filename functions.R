# =============================================================================
# boolean_helpers.R
#
# Utility functions for Boolean rule inference, attractor analysis,
# and GRN modeling. All functions use camelCase naming convention.
# =============================================================================

# =============================================================================
# timeIt
#
# =============================================================================

timeIt <- function(expr, name = "task") {
  start <- Sys.time()
  result <- force(expr)
  end   <- Sys.time()
  dur   <- difftime(end, start, units = "secs")
  cat(sprintf("%s finished in %d min %.1f s\n",
              name, as.integer(dur) %/% 60, as.numeric(dur) %% 60))
  invisible(result)
}
# -----------------------------------------------------------------------------
# Cell Type Menu
# -----------------------------------------------------------------------------

showCellTypeMenu <- function(cellTypes) {
  # Prompt the user only in an interactive R session
  if (is.null(cellTypes) || length(cellTypes) == 0) {
    stop("There are no valid cell types.")
  }
  if (interactive()) {
    cat("\014")
    cat("\n")
    selectionIndex <- menu(cellTypes, title = "Select a cell type:")
    
    if (selectionIndex == 0) {
      stop("No selection made. Exiting.")
    }
    cellType <- cellTypes[selectionIndex]
    cat("\014")
    message("You chose: ", cellType)
    message(" ")
  } else {
    stop("This script must be run interactively to choose a cell type.")
  }
  cellType
}

# -----------------------------------------------------------------------------
# Pseudotime Trajectory Menu
# -----------------------------------------------------------------------------

showTrajectoryMenu <- function(trajectories) {
  # Prompt the user only in an interactive R session
  if (is.null(trajectories) || length(trajectories) == 0) {
    stop("There are no valid pseudotime trajectories for this cell type.")
  }
  
  if (interactive()) {
    cat("\014")
    cat("\n")
    selectionIndex <- menu(trajectories, title = "Select a pseudotime trajectory:")
    
    if (selectionIndex == 0) {
      stop("No selection made. Exiting.")
    }
    traj <- trajectories[selectionIndex]
    cat("\014")
    message("You chose: ", traj)
    message(" ")
  } else {
    stop("This script must be run interactively to choose a cell type.")
  }
  traj
}

# -----------------------------------------------------------------------------
# Compute Input-Output Pairs for a Given Gene
# -----------------------------------------------------------------------------
makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k = 10) {
  ioList <- list()
  for (t in seq_len(length(cellOrder) - k)) {
    inputVec <- matBin[regulators, cellOrder[t], drop = FALSE]
    outputVal <- matBin[targetGene, cellOrder[t + k]]
    rowDf <- data.frame(t(inputVec), output = outputVal, stringsAsFactors = FALSE, check.names = FALSE)
    ioList[[t]] <- rowDf
  }
  if (length(ioList) == 0) return(NULL)
  do.call(rbind, ioList)
}

# -----------------------------------------------------------------------------
# Find the Best Boolean Rules from Input-Output Pairs
# -----------------------------------------------------------------------------
findBestBooleanRules <- function(ioDf) {
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  allInputStates <- expand.grid(rep(list(c(0,1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  allFnIdx <- seq_len(2^nStates) - 1L
  bestScore <- -1
  bestFns <- list()

  for (fnIdx in allFnIdx) {
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    correct <- 0
    for (i in seq_len(nrow(ioDf))) {
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == as.numeric(ioDf[i, inputCols]))))
      if (outPattern[matchIdx] == ioDf$output[i]) correct <- correct + 1
    }
    score <- correct / nrow(ioDf)
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
    } else if (abs(score - bestScore) < 1e-9) {
      bestFns[[length(bestFns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
    }
  }
  list(score = bestScore, bestFns = bestFns)
}

# -----------------------------------------------------------------------------
# Combine Best Boolean Rules Using Logical OR
# -----------------------------------------------------------------------------
combineBooleanFunctionsByOr <- function(fnList) {
  if (length(fnList) == 1) return(fnList[[1]]$outPattern)
  mat <- sapply(fnList, `[[`, "outPattern")
  as.integer(rowSums(mat) > 0)
}

# -----------------------------------------------------------------------------
# Generate BoolNet Rule String
# -----------------------------------------------------------------------------
makeBoolNetRule <- function(geneName, regulators, outPattern) {
  if (length(regulators) == 0)
    return(paste0(geneName, ", ", ifelse(outPattern[1] == 1, geneName, "! ", geneName)))

  R <- length(regulators)
  if (length(outPattern) != 2^R) stop("Mismatch: outPattern length != 2^numRegulators")

  clauseList <- c()
  for (i in 0:(2^R - 1)) {
    if (outPattern[i + 1] == 1) {
      bits <- as.integer(intToBits(i))[1:R]
      andTerms <- ifelse(bits == 1, regulators, paste0("! ", regulators))
      clauseList <- c(clauseList, paste0("(", paste(andTerms, collapse = " & "), ")"))
    }
  }
  if (length(clauseList) == 0) return(paste0(geneName, ", ! ", geneName))
  paste0(geneName, ", ", paste(clauseList, collapse = " | "))
}

# -----------------------------------------------------------------------------
# Gene Name Sanitization Helpers
# -----------------------------------------------------------------------------
sanitizeGeneName <- function(name) {
  name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  if (grepl("^[0-9]", name2)) name2 <- paste0("X", name2)
  name2
}

generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(OriginalName = geneNames, SanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Binary Vector Decoding from Big Integer
# -----------------------------------------------------------------------------
decodeBigIntegerState <- function(encodedVal, nGenes) {
  valBig <- as.bigz(encodedVal)
  bitVec <- numeric(nGenes)
  for (i in seq_len(nGenes)) {
    bitVec[i] <- as.integer(mod.bigz(valBig, as.bigz(2)))
    valBig <- valBig %/% as.bigz(2)
    if (valBig == 0) break
  }
  bitVec
}

# -----------------------------------------------------------------------------
# Shannon Entropy and Attractor Entropy Computation
# -----------------------------------------------------------------------------

shannonEntropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log2(p))
}

# -----------------------------------------------------------------------------

computeAttractorEntropy <- function(boolnet, attractors, nSamplesState = 5, nPerturb = 10) {
  geneNames <- attractors$stateInfo$genes
  nGenes <- length(geneNames)
  results <- data.frame(AttractorIndex = seq_along(attractors$attractors), Entropy = numeric(length(attractors$attractors)))
  
  for (i in seq_along(attractors$attractors)) {
    att <- attractors$attractors[[i]]
    allStates <- att$involvedStates
    if (length(allStates) == 0) next
    sel <- if (length(allStates) > nSamplesState) sample(allStates, nSamplesState) else allStates
    
    entropies <- sapply(sel, function(stateVal) {
      decoded <- decodeBigIntegerState(stateVal, nGenes)
      names(decoded) <- geneNames
      final <- replicate(nPerturb, {
        pert <- perturbNetwork(boolnet, method = "shuffle", perturb = "functions")
        getAttractors(pert, method = "chosen", startStates = list(decoded), returnTable = TRUE, type = "synchronous")$attractors[[1]]$involvedStates[[1]]
      })
      shannonEntropy(table(as.character(final)) / length(final))
    })
    results$Entropy[i] <- mean(entropies, na.rm = TRUE)
  }
  results$ScaledEntropy <- results$Entropy / log2(nrow(results))
  results$Stability <- 1 - results$ScaledEntropy
  results
}

# -----------------------------------------------------------------------------
# Function to test parameter convergence
# -----------------------------------------------------------------------------

testParameterConvergence <- function(boolnet, attractors, maxSamples = 50, maxPerturb = 50, 
                                     tolerance = 0.01, nReps = 5) {
  
  results <- data.frame(
    nSamplesState = integer(),
    nPerturb = integer(),
    MeanEntropy = numeric(),
    StdEntropy = numeric(),
    ComputeTime = numeric()
  )
  
  # Test different parameter combinations
  sampleValues <- c(5, 10, 15, 20, 25, 30, 40, 50)
  perturbValues <- c(5, 10, 15, 20, 25, 30, 40, 50)
  
  for (nSamp in sampleValues) {
    for (nPert in perturbValues) {
      if (nSamp > maxSamples || nPert > maxPerturb) next
      
      cat(sprintf("Testing nSamples=%d, nPerturb=%d\n", nSamp, nPert))
      
      # Run multiple replicates to assess stability
      entropies <- numeric(nReps)
      times <- numeric(nReps)
      
      for (rep in 1:nReps) {
        start_time <- Sys.time()
        
        tryCatch({
          entropy_result <- computeAttractorEntropy(boolnet, attractors, 
                                                    nSamplesState = nSamp, 
                                                    nPerturb = nPert)
          entropies[rep] <- mean(entropy_result$Entropy, na.rm = TRUE)
          times[rep] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        }, error = function(e) {
          entropies[rep] <- NA
          times[rep] <- NA
        })
      }
      
      # Store results
      results <- rbind(results, data.frame(
        nSamplesState = nSamp,
        nPerturb = nPert,
        MeanEntropy = mean(entropies, na.rm = TRUE),
        StdEntropy = sd(entropies, na.rm = TRUE),
        ComputeTime = mean(times, na.rm = TRUE)
      ))
    }
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Function to find optimal parameters based on convergence
# -----------------------------------------------------------------------------

findOptimalParameters <- function(convergence_results, time_weight = 0.3, stability_weight = 0.7) {
  
  # Normalize metrics (0-1 scale)
  convergence_results$NormTime <- 1 - (convergence_results$ComputeTime / max(convergence_results$ComputeTime, na.rm = TRUE))
  convergence_results$NormStability <- 1 - (convergence_results$StdEntropy / max(convergence_results$StdEntropy, na.rm = TRUE))
  
  # Composite score (higher is better)
  convergence_results$CompositeScore <- (time_weight * convergence_results$NormTime) + 
    (stability_weight * convergence_results$NormStability)
  
  # Find optimal parameters
  optimal_idx <- which.max(convergence_results$CompositeScore)
  optimal_params <- convergence_results[optimal_idx, ]
  
  return(list(
    optimal = optimal_params,
    all_results = convergence_results[order(-convergence_results$CompositeScore), ]
  ))
}

# -----------------------------------------------------------------------------
# Quick parameter test (much faster than full convergence analysis)
# -----------------------------------------------------------------------------

quickParameterTest <- function(boolnet, attractors) {
  
  cat("Running quick parameter optimization test...\n")
  cat("This will take about 2-3 hours but could save you days of compute time.\n\n")
  
  # Test just a few key combinations
  test_params <- data.frame(
    nSamplesState = c(10, 15, 20, 25),
    nPerturb = c(15, 20, 25, 25)
  )
  
  results <- data.frame(
    nSamplesState = integer(),
    nPerturb = integer(),
    MeanEntropy = numeric(),
    ComputeTime = numeric(),
    EstimatedFullTime = numeric()
  )
  
  for (i in 1:nrow(test_params)) {
    nSamp <- test_params$nSamplesState[i]
    nPert <- test_params$nPerturb[i]
    
    cat(sprintf("Testing (%d, %d)... ", nSamp, nPert))
    
    start_time <- Sys.time()
    
    tryCatch({
      entropy_result <- computeAttractorEntropy(boolnet, attractors, 
                                                nSamplesState = nSamp, 
                                                nPerturb = nPert)
      compute_time <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
      mean_entropy <- mean(entropy_result$Entropy, na.rm = TRUE)
      
      cat(sprintf("%.2f hours\n", compute_time))
      
      results <- rbind(results, data.frame(
        nSamplesState = nSamp,
        nPerturb = nPert,
        MeanEntropy = mean_entropy,
        ComputeTime = compute_time,
        EstimatedFullTime = compute_time  # This IS the full time for these params
      ))
      
    }, error = function(e) {
      cat("FAILED\n")
    })
  }
  
  # Add recommendations
  results$SpeedupVs25_25 <- 4.0 / results$ComputeTime  # Assuming 4 hours for (25,25)
  results$Recommended <- FALSE
  
  # Mark best options (good speedup with reasonable accuracy)
  if (nrow(results) > 0) {
    # Find parameters with >1.5x speedup
    good_speedup <- results$SpeedupVs25_25 >= 1.5
    if (any(good_speedup)) {
      best_idx <- which(good_speedup)[which.max(results$SpeedupVs25_25[good_speedup])]
      results$Recommended[best_idx] <- TRUE
    }
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Hybrid approach: Combine attractor and pseudotime analyses
# Add this to Script 13 for comprehensive target evaluation
# -----------------------------------------------------------------------------

evaluatePseudotimeShift <- function(cds, targetGenes, perturbType = "knockdown") {
  
  # Get baseline pseudotime distribution
  baseline_pt <- pseudotime(cds)
  baseline_mean <- mean(baseline_pt, na.rm = TRUE)
  
  # Simulate perturbation effects on pseudotime
  # This is a simplified simulation - in reality you'd need more sophisticated modeling
  
  results <- data.frame(
    Gene = targetGenes,
    BaselinePseudotime = baseline_mean,
    PredictedPseudotimeShift = numeric(length(targetGenes)),
    PseudotimeScore = numeric(length(targetGenes))
  )
  
  for (i in seq_along(targetGenes)) {
    gene <- targetGenes[i]
    
    # Simple correlation-based prediction of pseudotime shift
    # In practice, you might use more sophisticated methods
    gene_expr <- assay(cds, "logcounts")[gene, ]
    pt_correlation <- cor(gene_expr, baseline_pt, use = "complete.obs")
    
    # Predict pseudotime shift based on perturbation direction
    if (perturbType == "knockdown") {
      predicted_shift <- -pt_correlation * sd(baseline_pt) * 0.5  # Heuristic
    } else {  # overexpression
      predicted_shift <- pt_correlation * sd(baseline_pt) * 0.5
    }
    
    results$PredictedPseudotimeShift[i] <- predicted_shift
    
    # Score: negative shifts (toward young) are better
    results$PseudotimeScore[i] <- -predicted_shift / sd(baseline_pt)
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Comprehensive target evaluation combining both approaches
# -----------------------------------------------------------------------------

comprehensiveTargetEvaluation <- function(attractorResults, cds, targetGenes) {
  
  # Get attractor-based scores
  attractor_scores <- attractorResults %>%
    arrange(AgingScore) %>%
    mutate(AttractorRank = row_number(),
           AttractorScore = (nrow(.) - AttractorRank + 1) / nrow(.))  # Normalize to 0-1
  
  # Get pseudotime-based scores
  pt_kd_scores <- evaluatePseudotimeShift(cds, targetGenes, "knockdown")
  pt_oe_scores <- evaluatePseudotimeShift(cds, targetGenes, "overexpression")
  
  # Combine results
  combined <- attractor_scores %>%
    left_join(pt_kd_scores %>% select(Gene, PseudotimeScore_KD = PseudotimeScore), by = "Gene") %>%
    left_join(pt_oe_scores %>% select(Gene, PseudotimeScore_OE = PseudotimeScore), by = "Gene") %>%
    mutate(
      # Choose best pseudotime score for each gene
      BestPseudotimeScore = pmax(PseudotimeScore_KD, PseudotimeScore_OE, na.rm = TRUE),
      
      # Composite score combining both approaches
      CompositeScore = 0.6 * AttractorScore + 0.4 * BestPseudotimeScore,
      
      # Agreement metric: do both approaches agree on this target?
      Agreement = case_when(
        AttractorScore > 0.7 & BestPseudotimeScore > 0.7 ~ "High Agreement",
        AttractorScore > 0.5 & BestPseudotimeScore > 0.5 ~ "Moderate Agreement", 
        (AttractorScore > 0.7) != (BestPseudotimeScore > 0.7) ~ "Disagreement",
        TRUE ~ "Low Agreement"
      )
    ) %>%
    arrange(desc(CompositeScore))
  
  return(combined)
}

# -----------------------------------------------------------------------------
# Target validation priority system
# -----------------------------------------------------------------------------

prioritizeTargetsForValidation <- function(comprehensive_results) {
  
  priorities <- comprehensive_results %>%
    mutate(
      Priority = case_when(
        Agreement == "High Agreement" & CompositeScore > 0.8 ~ "Tier 1: High Priority",
        Agreement == "High Agreement" & CompositeScore > 0.6 ~ "Tier 2: Medium Priority",
        Agreement == "Moderate Agreement" & CompositeScore > 0.7 ~ "Tier 2: Medium Priority",
        Agreement == "Disagreement" & AttractorScore > 0.8 ~ "Tier 3: Attractor-Specific",
        Agreement == "Disagreement" & BestPseudotimeScore > 0.8 ~ "Tier 3: Trajectory-Specific",
        TRUE ~ "Tier 4: Low Priority"
      )
    ) %>%
    arrange(match(Priority, c("Tier 1: High Priority", "Tier 2: Medium Priority", 
                              "Tier 3: Attractor-Specific", "Tier 3: Trajectory-Specific", 
                              "Tier 4: Low Priority")), desc(CompositeScore))
  
  return(priorities)
}

# -----------------------------------------------------------------------------
# Optimized Boolean rule inference with multiple strategies
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Strategy 1: Early termination for perfect rules
# -----------------------------------------------------------------------------

findBestBooleanRulesOptimized <- function(ioDf, early_termination = TRUE, max_functions = 10000) {
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  allInputStates <- expand.grid(rep(list(c(0,1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  allFnIdx <- seq_len(2^nStates) - 1L
  
  # Early termination if too many functions to test
  if (length(allFnIdx) > max_functions) {
    warning(sprintf("Too many functions (%d) for exhaustive search. Using sampling approach.", 
                    length(allFnIdx)))
    return(findBooleanRulesSampling(ioDf, max_functions))
  }
  
  bestScore <- -1
  bestFns <- list()
  perfectFound <- FALSE
  
  for (fnIdx in allFnIdx) {
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    correct <- 0
    
    for (i in seq_len(nrow(ioDf))) {
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == as.numeric(ioDf[i, inputCols]))))
      if (outPattern[matchIdx] == ioDf$output[i]) correct <- correct + 1
    }
    
    score <- correct / nrow(ioDf)
    
    if (score > bestScore) {
      bestScore <- score
      bestFns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
      
      # Early termination for perfect score
      if (early_termination && score >= 0.999) {
        perfectFound <- TRUE
        break
      }
    } else if (abs(score - bestScore) < 1e-9) {
      bestFns[[length(bestFns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
    }
  }
  
  return(list(score = bestScore, bestFns = bestFns, perfectFound = perfectFound))
}

# -----------------------------------------------------------------------------
# Strategy 2: Intelligent sampling for high-regulator cases
# -----------------------------------------------------------------------------

findBooleanRulesSampling <- function(ioDf, max_samples = 10000, method = "smart") {
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  allInputStates <- expand.grid(rep(list(c(0,1)), R))
  colnames(allInputStates) <- inputCols
  nStates <- 2^R
  totalFunctions <- 2^nStates
  
  if (method == "smart") {
    # Bias sampling toward simpler functions (fewer 1s in truth table)
    candidates <- generateSmartCandidates(ioDf, allInputStates, max_samples)
  } else {
    # Pure random sampling
    candidates <- sample(0:(totalFunctions-1), min(max_samples, totalFunctions))
  }
  
  bestScore <- -1
  bestFns <- list()
  
  for (fnIdx in candidates) {
    outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
    correct <- 0
    
    for (i in seq_len(nrow(ioDf))) {
      matchIdx <- which(apply(allInputStates, 1, function(x) all(x == as.numeric(ioDf[i, inputCols]))))
      if (outPattern[matchIdx] == ioDf$output[i]) correct <- correct + 1
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

# -----------------------------------------------------------------------------
# Strategy 3: Generate biologically-plausible candidates
# -----------------------------------------------------------------------------

generateSmartCandidates <- function(ioDf, allInputStates, max_samples) {
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  nStates <- 2^R
  
  candidates <- c()
  
  # 1. Simple functions (AND, OR, single inputs)
  simple_functions <- c(
    0,                    # Always 0
    2^nStates - 1,       # Always 1
    as.integer(paste0(rep("1", nStates), collapse = ""), base = 2)  # All combinations
  )
  
  # 2. Single regulator functions
  for (reg in seq_along(inputCols)) {
    pattern <- rep(0, nStates)
    for (state in 1:nStates) {
      bits <- as.integer(intToBits(state-1))[1:R]
      if (bits[reg] == 1) pattern[state] <- 1
    }
    fnIdx <- sum(pattern * 2^(0:(nStates-1)))
    candidates <- c(candidates, fnIdx)
    
    # Negation
    pattern_neg <- 1 - pattern
    fnIdx_neg <- sum(pattern_neg * 2^(0:(nStates-1)))
    candidates <- c(candidates, fnIdx_neg)
  }
  
  # 3. Two-regulator combinations (AND, OR, XOR)
  if (R >= 2) {
    for (i in 1:(R-1)) {
      for (j in (i+1):R) {
        # AND
        pattern_and <- rep(0, nStates)
        for (state in 1:nStates) {
          bits <- as.integer(intToBits(state-1))[1:R]
          if (bits[i] == 1 && bits[j] == 1) pattern_and[state] <- 1
        }
        candidates <- c(candidates, sum(pattern_and * 2^(0:(nStates-1))))
        
        # OR
        pattern_or <- rep(0, nStates)
        for (state in 1:nStates) {
          bits <- as.integer(intToBits(state-1))[1:R]
          if (bits[i] == 1 || bits[j] == 1) pattern_or[state] <- 1
        }
        candidates <- c(candidates, sum(pattern_or * 2^(0:(nStates-1))))
      }
    }
  }
  
  # 4. Random sampling for remaining slots
  remaining <- max_samples - length(unique(candidates))
  if (remaining > 0) {
    random_candidates <- sample(0:(2^nStates - 1), remaining)
    candidates <- c(candidates, random_candidates)
  }
  
  return(unique(candidates)[1:min(max_samples, length(unique(candidates)))])
}

# -----------------------------------------------------------------------------
# Strategy 4: Parallel optimization
# -----------------------------------------------------------------------------

findBooleanRulesParallel <- function(ioDf, cores = NULL) {
  if (is.null(cores)) cores <- min(4, parallel::detectCores() - 1)
  
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  nStates <- 2^R
  allFnIdx <- seq_len(2^nStates) - 1L
  
  # Split function indices across cores
  chunks <- split(allFnIdx, cut(seq_along(allFnIdx), cores))
  
  # Parallel computation
  chunk_results <- foreach(chunk = chunks, .combine = 'c', .packages = c()) %dopar% {
    chunk_best <- list(score = -1, fns = list())
    allInputStates <- expand.grid(rep(list(c(0,1)), R))
    colnames(allInputStates) <- inputCols
    
    for (fnIdx in chunk) {
      outPattern <- as.integer(intToBits(fnIdx)[1:nStates])
      correct <- 0
      
      for (i in seq_len(nrow(ioDf))) {
        matchIdx <- which(apply(allInputStates, 1, function(x) all(x == as.numeric(ioDf[i, inputCols]))))
        if (outPattern[matchIdx] == ioDf$output[i]) correct <- correct + 1
      }
      
      score <- correct / nrow(ioDf)
      
      if (score > chunk_best$score) {
        chunk_best$score <- score
        chunk_best$fns <- list(list(patternIdx = fnIdx, outPattern = outPattern, score = score))
      } else if (abs(score - chunk_best$score) < 1e-9) {
        chunk_best$fns[[length(chunk_best$fns) + 1]] <- list(patternIdx = fnIdx, outPattern = outPattern, score = score)
      }
    }
    
    list(chunk_best)
  }
  
  # Combine results from all chunks
  bestScore <- max(sapply(chunk_results, function(x) x$score))
  bestFns <- do.call(c, lapply(chunk_results, function(x) {
    if (abs(x$score - bestScore) < 1e-9) x$fns else list()
  }))
  
  return(list(score = bestScore, bestFns = bestFns))
}

# -----------------------------------------------------------------------------
# Strategy 5: Adaptive approach that chooses method based on problem size
# -----------------------------------------------------------------------------

findBestBooleanRulesAdaptive <- function(ioDf, time_limit_minutes = 30) {
  inputCols <- setdiff(colnames(ioDf), "output")
  R <- length(inputCols)
  totalFunctions <- 2^(2^R)
  
  start_time <- Sys.time()
  
  if (R <= 2) {
    # Small problem: use exact method
    return(findBestBooleanRulesOptimized(ioDf, early_termination = TRUE))
    
  } else if (R == 3) {
    # Medium problem: try exact with early termination, fall back to parallel
    result <- findBestBooleanRulesOptimized(ioDf, early_termination = TRUE)
    if (result$perfectFound || difftime(Sys.time(), start_time, units = "mins") < time_limit_minutes/2) {
      return(result)
    } else {
      return(findBooleanRulesParallel(ioDf))
    }
    
  } else if (R <= 4) {
    # Large problem: use smart sampling
    return(findBooleanRulesSampling(ioDf, max_samples = 50000, method = "smart"))
    
  } else {
    # Very large problem: use aggressive sampling
    return(findBooleanRulesSampling(ioDf, max_samples = 10000, method = "smart"))
  }
}
