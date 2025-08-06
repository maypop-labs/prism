# =============================================================================
# attractorManager.R
# Purpose: Functions for attractor analysis, entropy computation, and target evaluation
# Dependencies: BoolNet, gmp, monocle3
# =============================================================================

# =============================================================================
# Gene Name Sanitization
# =============================================================================

#' Sanitize gene names for BoolNet compatibility
#'
#' Converts gene names to valid identifiers by replacing special characters
#' with underscores and ensuring names don't start with numbers.
#'
#' @param name Character string representing a gene name
#' @return Sanitized gene name suitable for BoolNet
#' @export
#' @examples
#' sanitizeGeneName("IL-6")     # Returns "IL_6"
#' sanitizeGeneName("3-PGDH")   # Returns "X3_PGDH"
sanitizeGeneName <- function(name) {
  # Remove leading/trailing whitespace and replace special chars with underscore
  name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  
  # Prepend 'X' if name starts with a number
  if (grepl("^[0-9]", name2)) {
    name2 <- paste0("X", name2)
  }
  
  return(name2)
}

#' Generate mapping between original and sanitized gene names
#'
#' Creates a data frame mapping original gene names to their sanitized versions
#' for later reference and conversion.
#'
#' @param geneNames Character vector of original gene names
#' @return Data frame with columns 'OriginalName' and 'SanitizedName'
#' @export
generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(
    originalName = geneNames,
    sanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Binary State Decoding
# =============================================================================

#' Decode big integer state to binary vector
#'
#' Converts a BoolNet attractor state (encoded as big integer) back to
#' a binary vector representing gene ON/OFF states.
#'
#' @param encodedVal Big integer encoded state from BoolNet
#' @param nGenes Number of genes in the network
#' @return Numeric vector of 0s and 1s representing gene states
#' @export
#' @importFrom gmp as.bigz mod.bigz
decodeBigIntegerState <- function(encodedVal, nGenes) {
  valBig <- as.bigz(encodedVal)
  bitVec <- numeric(nGenes)
  
  for (i in seq_len(nGenes)) {
    bitVec[i] <- as.integer(mod.bigz(valBig, as.bigz(2)))
    valBig <- valBig %/% as.bigz(2)
    
    # Early termination if no more bits to process
    if (valBig == 0) break
  }
  
  return(bitVec)
}

# =============================================================================
# Entropy Calculations
# =============================================================================

#' Calculate Shannon entropy from probability vector
#'
#' Computes Shannon entropy H = -sum(p * log2(p)) for a probability distribution.
#' Handles zero probabilities gracefully.
#'
#' @param p Numeric vector of probabilities (should sum to 1)
#' @return Shannon entropy value in bits
#' @export
#' @examples
#' shannonEntropy(c(0.5, 0.5))    # Returns 1.0 (maximum for 2 states)
#' shannonEntropy(c(1.0, 0.0))    # Returns 0.0 (no uncertainty)
shannonEntropy <- function(p) {
  # Remove zero probabilities to avoid log(0)
  p <- p[p > 0]
  
  if (length(p) == 0) return(0)
  
  return(-sum(p * log2(p)))
}

# =============================================================================
# Attractor Entropy Analysis
# =============================================================================

#' Compute entropy and stability metrics for Boolean network attractors
#'
#' Analyzes the stability of each attractor by perturbing the network and
#' measuring how often the system returns to the same attractor state.
#' Higher entropy indicates less stable attractors.
#'
#' @param boolnet BoolNet network object
#' @param attractors Attractor object from getAttractors()
#' @param nSamplesState Number of states to sample per attractor (default: 5)
#' @param nPerturb Number of perturbation trials per state (default: 10)
#' @return Data frame with columns: AttractorIndex, Entropy, ScaledEntropy, Stability
#' @export
#' @importFrom BoolNet perturbNetwork getAttractors
computeAttractorEntropy <- function(boolnet, attractors, nSamplesState = 5, nPerturb = 10) {
  geneNames <- attractors$stateInfo$genes
  nGenes <- length(geneNames)
  nAttractors <- length(attractors$attractors)
  
  # Initialize results data frame
  results <- data.frame(
    AttractorIndex = seq_along(attractors$attractors),
    Entropy = numeric(nAttractors)
  )
  
  # Analyze each attractor
  for (i in seq_along(attractors$attractors)) {
    att <- attractors$attractors[[i]]
    allStates <- att$involvedStates
    
    # Skip empty attractors
    if (length(allStates) == 0) {
      results$Entropy[i] <- NA
      next
    }
    
    # Sample states from this attractor
    statesToTest <- if (length(allStates) > nSamplesState) {
      sample(allStates, nSamplesState)
    } else {
      allStates
    }
    
    # Calculate entropy for each sampled state
    stateEntropies <- sapply(statesToTest, function(stateVal) {
      # Decode state to binary vector
      decoded <- decodeBigIntegerState(stateVal, nGenes)
      names(decoded) <- geneNames
      
      # Perform perturbation experiments
      finalStates <- replicate(nPerturb, {
        tryCatch({
          # Perturb network and find new attractor
          perturbedNet <- perturbNetwork(boolnet, method = "shuffle", perturb = "functions")
          newAttractors <- getAttractors(perturbedNet, 
                                         method = "chosen", 
                                         startStates = list(decoded), 
                                         returnTable = TRUE, 
                                         type = "synchronous")
          
          # Return final state
          newAttractors$attractors[[1]]$involvedStates[[1]]
        }, error = function(e) {
          # Return original state if perturbation fails
          stateVal
        })
      })
      
      # Calculate entropy of final state distribution
      stateCounts <- table(as.character(finalStates))
      probabilities <- stateCounts / length(finalStates)
      return(shannonEntropy(probabilities))
    })
    
    # Average entropy across all sampled states
    results$Entropy[i] <- mean(stateEntropies, na.rm = TRUE)
  }
  
  # Add derived metrics
  maxEntropy <- log2(nrow(results))
  results$ScaledEntropy <- results$Entropy / maxEntropy
  results$Stability <- 1 - results$ScaledEntropy
  
  return(results)
}

# =============================================================================
# Parameter Optimization
# =============================================================================

#' Test different parameter combinations for entropy computation
#'
#' Evaluates various nSamplesState and nPerturb parameter combinations to find
#' optimal balance between computational cost and result stability.
#'
#' @param boolnet BoolNet network object
#' @param attractors Attractor object from getAttractors()
#' @param maxSamples Maximum number of samples to test (default: 50)
#' @param maxPerturb Maximum number of perturbations to test (default: 50)
#' @param tolerance Convergence tolerance (default: 0.01)
#' @param nReps Number of replicates per parameter combination (default: 5)
#' @return Data frame with parameter combinations and their performance metrics
#' @export
testParameterConvergence <- function(boolnet, attractors, maxSamples = 50, maxPerturb = 50, 
                                     tolerance = 0.01, nReps = 5) {
  
  results <- data.frame(
    nSamplesState = integer(),
    nPerturb = integer(),
    MeanEntropy = numeric(),
    StdEntropy = numeric(),
    ComputeTime = numeric()
  )
  
  # Define parameter ranges to test
  sampleValues <- c(5, 10, 15, 20, 25, 30, 40, 50)
  perturbValues <- c(5, 10, 15, 20, 25, 30, 40, 50)
  
  message("Testing parameter convergence...")
  
  for (nSamp in sampleValues) {
    for (nPert in perturbValues) {
      # Skip if parameters exceed limits
      if (nSamp > maxSamples || nPert > maxPerturb) next
      
      message(sprintf("Testing nSamples=%d, nPerturb=%d", nSamp, nPert))
      
      # Run multiple replicates
      entropies <- numeric(nReps)
      times <- numeric(nReps)
      
      for (rep in 1:nReps) {
        startTime <- Sys.time()
        
        tryCatch({
          entropyResult <- computeAttractorEntropy(boolnet, attractors, 
                                                   nSamplesState = nSamp, 
                                                   nPerturb = nPert)
          entropies[rep] <- mean(entropyResult$Entropy, na.rm = TRUE)
          times[rep] <- as.numeric(difftime(Sys.time(), startTime, units = "secs"))
        }, error = function(e) {
          entropies[rep] <- NA
          times[rep] <- NA
          warning(sprintf("Failed for parameters (%d, %d): %s", nSamp, nPert, e$message))
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

#' Find optimal parameters based on convergence analysis
#'
#' Analyzes parameter convergence results to recommend optimal settings
#' balancing computational efficiency and result stability.
#'
#' @param convergenceResults Data frame from testParameterConvergence()
#' @param timeWeight Weight for computational time in scoring (default: 0.3)
#' @param stabilityWeight Weight for result stability in scoring (default: 0.7)
#' @return List with optimal parameters and ranked results
#' @export
findOptimalParameters <- function(convergenceResults, timeWeight = 0.3, stabilityWeight = 0.7) {
  
  # Remove rows with missing data
  validResults <- convergenceResults[complete.cases(convergenceResults), ]
  
  if (nrow(validResults) == 0) {
    stop("No valid convergence results found")
  }
  
  # Normalize metrics to 0-1 scale (higher is better)
  validResults$NormTime <- 1 - (validResults$ComputeTime / max(validResults$ComputeTime, na.rm = TRUE))
  validResults$NormStability <- 1 - (validResults$StdEntropy / max(validResults$StdEntropy, na.rm = TRUE))
  
  # Calculate composite score
  validResults$CompositeScore <- (timeWeight * validResults$NormTime) + 
    (stabilityWeight * validResults$NormStability)
  
  # Find optimal parameters
  optimalIdx <- which.max(validResults$CompositeScore)
  optimalParams <- validResults[optimalIdx, ]
  
  # Sort all results by score
  rankedResults <- validResults[order(-validResults$CompositeScore), ]
  
  return(list(
    optimal = optimalParams,
    all_results = rankedResults
  ))
}

#' Quick parameter optimization test
#'
#' Performs a faster parameter test using predefined combinations to give
#' quick recommendations without extensive computation.
#'
#' @param boolnet BoolNet network object
#' @param attractors Attractor object from getAttractors()
#' @return Data frame with parameter recommendations and estimated compute times
#' @export
quickParameterTest <- function(boolnet, attractors) {
  
  message("Running quick parameter optimization test...")
  message("This will take 2-3 hours but could save days of compute time.\n")
  
  # Test predefined parameter combinations
  testParams <- data.frame(
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
  
  for (i in 1:nrow(testParams)) {
    nSamp <- testParams$nSamplesState[i]
    nPert <- testParams$nPerturb[i]
    
    message(sprintf("Testing (%d, %d)... ", nSamp, nPert), appendLF = FALSE)
    
    startTime <- Sys.time()
    
    tryCatch({
      entropyResult <- computeAttractorEntropy(boolnet, attractors, 
                                               nSamplesState = nSamp, 
                                               nPerturb = nPert)
      computeTime <- as.numeric(difftime(Sys.time(), startTime, units = "hours"))
      meanEntropy <- mean(entropyResult$Entropy, na.rm = TRUE)
      
      message(sprintf("%.2f hours", computeTime))
      
      results <- rbind(results, data.frame(
        nSamplesState = nSamp,
        nPerturb = nPert,
        MeanEntropy = meanEntropy,
        ComputeTime = computeTime,
        EstimatedFullTime = computeTime
      ))
      
    }, error = function(e) {
      message("FAILED")
      warning(sprintf("Parameter test failed for (%d, %d): %s", nSamp, nPert, e$message))
    })
  }
  
  # Add performance metrics
  if (nrow(results) > 0) {
    # Estimate speedup compared to most expensive option
    maxTime <- max(results$ComputeTime, na.rm = TRUE)
    results$SpeedupVs25_25 <- maxTime / results$ComputeTime
    results$Recommended <- FALSE
    
    # Mark best option (good speedup with reasonable performance)
    goodSpeedup <- results$SpeedupVs25_25 >= 1.5
    if (any(goodSpeedup)) {
      bestIdx <- which(goodSpeedup)[which.max(results$SpeedupVs25_25[goodSpeedup])]
      results$Recommended[bestIdx] <- TRUE
    }
  }
  
  return(results)
}

# =============================================================================
# Comprehensive Target Evaluation
# =============================================================================

#' Evaluate pseudotime shift predictions for target genes
#'
#' Predicts how perturbing target genes might affect cellular pseudotime
#' progression based on gene-pseudotime correlations.
#'
#' @param cds Monocle3 CellDataSet object
#' @param targetGenes Character vector of gene names to evaluate
#' @param perturbType Type of perturbation ("knockdown" or "overexpression")
#' @return Data frame with predicted pseudotime shifts and scores
#' @export
evaluatePseudotimeShift <- function(cds, targetGenes, perturbType = "knockdown") {
  
  # Get baseline pseudotime distribution
  baselinePt <- pseudotime(cds)
  baselineMean <- mean(baselinePt, na.rm = TRUE)
  
  results <- data.frame(
    Gene = targetGenes,
    BaselinePseudotime = baselineMean,
    PredictedPseudotimeShift = numeric(length(targetGenes)),
    PseudotimeScore = numeric(length(targetGenes))
  )
  
  for (i in seq_along(targetGenes)) {
    gene <- targetGenes[i]
    
    # Check if gene exists in dataset
    if (!gene %in% rownames(cds)) {
      warning(sprintf("Gene %s not found in dataset", gene))
      results$PredictedPseudotimeShift[i] <- NA
      results$PseudotimeScore[i] <- NA
      next
    }
    
    # Calculate gene-pseudotime correlation
    geneExpr <- assay(cds, "logcounts")[gene, ]
    ptCorrelation <- cor(geneExpr, baselinePt, use = "complete.obs")
    
    # Predict pseudotime shift based on perturbation direction
    if (perturbType == "knockdown") {
      predictedShift <- -ptCorrelation * sd(baselinePt, na.rm = TRUE) * 0.5
    } else {  # overexpression
      predictedShift <- ptCorrelation * sd(baselinePt, na.rm = TRUE) * 0.5
    }
    
    results$PredictedPseudotimeShift[i] <- predictedShift
    
    # Score: negative shifts (toward young) are better
    results$PseudotimeScore[i] <- -predictedShift / sd(baselinePt, na.rm = TRUE)
  }
  
  return(results)
}

#' Comprehensive target evaluation combining attractor and pseudotime analyses
#'
#' Integrates multiple evaluation approaches to rank potential therapeutic targets
#' based on both Boolean network attractor analysis and pseudotime predictions.
#'
#' @param attractorResults Data frame with attractor-based target scores
#' @param cds Monocle3 CellDataSet object
#' @param targetGenes Character vector of genes to evaluate
#' @return Data frame with comprehensive target rankings
#' @export
#' @importFrom dplyr arrange mutate left_join
comprehensiveTargetEvaluation <- function(attractorResults, cds, targetGenes) {
  
  # Get attractor-based scores
  attractorScores <- attractorResults %>%
    arrange(AgingScore) %>%
    mutate(AttractorRank = row_number(),
           AttractorScore = (nrow(.) - AttractorRank + 1) / nrow(.))
  
  # Get pseudotime-based scores
  ptKdScores <- evaluatePseudotimeShift(cds, targetGenes, "knockdown")
  ptOeScores <- evaluatePseudotimeShift(cds, targetGenes, "overexpression")
  
  # Combine results
  combined <- attractorScores %>%
    left_join(ptKdScores %>% select(Gene, PseudotimeScore_KD = PseudotimeScore), by = "Gene") %>%
    left_join(ptOeScores %>% select(Gene, PseudotimeScore_OE = PseudotimeScore), by = "Gene") %>%
    mutate(
      # Choose best pseudotime score for each gene
      BestPseudotimeScore = pmax(PseudotimeScore_KD, PseudotimeScore_OE, na.rm = TRUE),
      
      # Composite score combining both approaches
      CompositeScore = 0.6 * AttractorScore + 0.4 * BestPseudotimeScore,
      
      # Agreement metric between approaches
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

#' Prioritize targets for experimental validation
#'
#' Creates a tiered prioritization system for target validation based on
#' agreement between different analytical approaches and composite scores.
#'
#' @param comprehensiveResults Data frame from comprehensiveTargetEvaluation()
#' @return Data frame with validation priority tiers
#' @export
#' @importFrom dplyr mutate arrange
prioritizeTargetsForValidation <- function(comprehensiveResults) {
  
  priorities <- comprehensiveResults %>%
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
                              "Tier 4: Low Priority")), 
            desc(CompositeScore))
  
  return(priorities)
}