# =============================================================================
# attractorManager.R
# Purpose: Functions for attractor analysis and binary state decoding
# Dependencies: BoolNet, gmp
# =============================================================================

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
# Attractor Entropy and Stability Analysis
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
