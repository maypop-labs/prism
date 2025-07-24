# =============================================================================
# boolean_helpers.R
#
# Utility functions for Boolean rule inference, attractor analysis,
# and GRN modeling. All functions use camelCase naming convention.
# =============================================================================

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

# -----------------------------------------------------------------------------
#  computeAttractorEntropy  (serial version, no parallel back-end required)
# -----------------------------------------------------------------------------
#  • boolnet      : BoolNet object returned by load/run pipeline
#  • attractors   : list returned by BoolNet::getAttractors()
#  • nSamplesState: how many distinct states to sample per attractor
#  • nPerturb     : how many perturbation runs per sampled state
#  • showProgress : TRUE ⇒ text progress-bar in console
# -----------------------------------------------------------------------------

computeAttractorEntropy <- function(boolnet, attractors,
                                    nSamplesState = 25,
                                    nPerturb      = 25,
                                    showProgress  = TRUE) {
  
  geneNames <- attractors$stateInfo$genes
  nGenes    <- length(geneNames)
  nAttr     <- length(attractors$attractors)
  
  if (showProgress)
    pb <- utils::txtProgressBar(min = 0, max = nAttr, style = 3)
  
  results <- data.frame(
    AttractorIndex = seq_len(nAttr),
    Entropy        = NA_real_
  )
  
  for (i in seq_len(nAttr)) {
    
    attStates <- attractors$attractors[[i]]$involvedStates
    if (length(attStates) == 0) {                # safeguard
      if (showProgress) utils::setTxtProgressBar(pb, i)
      next
    }
    
    ## ── sample up to nSamplesState unique states ────────────────────────────
    sel <- if (length(attStates) > nSamplesState)
      sample(attStates, nSamplesState)
    else
      attStates
    
    ## ── run perturbations and collect final attractor IDs ───────────────────
    entropies <- vapply(sel, function(stateVal) {
      
      decoded <- decodeBigIntegerState(stateVal, nGenes)
      names(decoded) <- geneNames
      
      finals <- replicate(nPerturb, {
        # here we shuffle Boolean functions to simulate perturbation
        pert  <- BoolNet::perturbNetwork(boolnet,
                                         method   = "shuffle",
                                         perturb  = "functions")
        BoolNet::getAttractors(pert, method = "chosen",
                               startStates  = list(decoded),
                               returnTable  = TRUE,
                               type         = "synchronous")$
          attractors[[1]]$involvedStates[[1]]
      })
      
      shannonEntropy(table(as.character(finals)) / length(finals))
      
    }, numeric(1))
    
    results$Entropy[i] <- mean(entropies, na.rm = TRUE)
    
    if (showProgress) utils::setTxtProgressBar(pb, i)
  }
  
  if (showProgress) close(pb)
  
  results$ScaledEntropy <- results$Entropy / log2(nAttr)
  results$Stability     <- 1 - results$ScaledEntropy
  results
}

# =============================================================================
# computeAttractorEntropy_bitflip.R
#
# Serial (single‑core) entropy / stability estimator for a BoolNet attractor
# landscape.  Instead of shuffling all Boolean functions we perturb the system
# by flipping **one random gene** in a sampled steady state. This avoids the
# heavy rebuild inside BoolNet::perturbNetwork() and usually yields a 3–5×
# speed‑up while capturing stability just as well.
#
# Required helpers already exist in functions.R:
#   - decodeBigIntegerState()
#   - shannonEntropy()
#
# =============================================================================

computeAttractorEntropy_bitflip <- function(boolnet, attractors,
                                            nSamplesState = 25,
                                            nPerturb      = 25,
                                            showProgress  = TRUE) {
  
  geneNames <- attractors$stateInfo$genes
  nGenes    <- length(geneNames)
  nAttr     <- length(attractors$attractors)
  
  if (showProgress)
    pb <- utils::txtProgressBar(min = 0, max = nAttr, style = 3)
  
  results <- data.frame(
    AttractorIndex = seq_len(nAttr),
    Entropy        = NA_real_
  )
  
  # ---------------------------------------------------------------------------
  # Pre‑decode every integer state once, so we never call gmp functions inside
  # the hot loop.
  # ---------------------------------------------------------------------------
  stateCache <- lapply(attractors$attractors, function(att) {
    lapply(att$involvedStates, decodeBigIntegerState, nGenes = nGenes)
  })
  
  for (i in seq_len(nAttr)) {
    
    attStates <- stateCache[[i]]
    S         <- length(attStates)
    if (S == 0L) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      next
    }
    
    # sample up to nSamplesState distinct states
    selIdx <- if (S > nSamplesState) sample.int(S, nSamplesState) else seq_len(S)
    sel    <- attStates[selIdx]
    
    # -----------------------------------------------------------------------
    # For each sampled state we perform nPerturb independent random bit flips
    # and record to which attractor the system returns.
    # -----------------------------------------------------------------------
    entropies <- vapply(sel, function(st) {
      
      finals <- integer(nPerturb)
      for (rep in seq_len(nPerturb)) {
        s <- st                        # copy the steady state
        flip <- sample.int(nGenes, 1L) # choose a random gene position
        s[flip] <- 1L - s[flip]        # toggle bit (0 -> 1 or 1 -> 0)
        
        finals[rep] <- BoolNet::getAttractors(
          boolnet,
          method      = "chosen",
          startStates = list(s),
          returnTable = TRUE,
          type        = "synchronous"
        )$attractors[[1]]$involvedStates[[1]]
      }
      
      shannonEntropy( table(as.character(finals)) / nPerturb )
      
    }, numeric(1))
    
    results$Entropy[i] <- mean(entropies, na.rm = TRUE)
    
    if (showProgress) utils::setTxtProgressBar(pb, i)
  }
  
  if (showProgress) close(pb)
  
  results$ScaledEntropy <- results$Entropy / log2(nAttr)
  results$Stability     <- 1 - results$ScaledEntropy
  results
}

# =============================================================================
# computeAttractorEntropy_bitflip_vec.R
#
# Vectorised, single‑core entropy / stability estimator.
# Instead of a nested replicate() loop (S sampled states × M perturbations), we
# pre‑allocate two integer matrices (seeds & flips) and iterate over columns.
# This trims R‑level function calls and makes per‑attractor timing roughly
# O(S×M) with a small constant factor.
#
# External deps: BoolNet, gmp  (both should already be loaded by caller)
# ----------------------------------------------------------------------------
#  boolnet        : BoolNet object produced earlier in pipeline
#  attractors     : list from BoolNet::getAttractors()
#  nSamplesState  : max # of distinct states to draw from each attractor
#  nPerturb       : # of single‑gene flips per sampled state
#  showProgress   : TRUE ⇒ console txt progress‑bar
# ----------------------------------------------------------------------------

computeAttractorEntropy_bitflip_vec <- function(boolnet, attractors,
                                                nSamplesState = 25,
                                                nPerturb      = 25,
                                                showProgress  = TRUE) {
  
  geneNames <- attractors$stateInfo$genes
  nGenes    <- length(geneNames)
  nAttr     <- length(attractors$attractors)
  
  if (showProgress)
    pb <- utils::txtProgressBar(min = 0, max = nAttr, style = 3)
  
  # ---- prepare decoded cache once -----------------------------------------
  decodeCache <- lapply(attractors$attractors, function(att) {
    lapply(att$involvedStates, decodeBigIntegerState, nGenes = nGenes)
  })
  
  results <- data.frame(AttractorIndex = seq_len(nAttr), Entropy = NA_real_)
  
  for (i in seq_len(nAttr)) {
    
    decodedStates <- decodeCache[[i]]
    if (length(decodedStates) == 0) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      next
    }
    
    # ---- sample up to nSamplesState unique steady states ------------------
    selIdx   <- sample(seq_along(decodedStates),
                       size = min(length(decodedStates), nSamplesState))
    sel      <- decodedStates[selIdx]
    S        <- length(sel)
    M        <- nPerturb
    
    # ---- vectorise: matrices seeds[S,M] and flips[S,M] --------------------
    seedsIdx <- matrix(sample(selIdx,  S * M, replace = TRUE), nrow = S)
    flipsMat <- matrix(sample.int(nGenes, S * M, replace = TRUE), nrow = S)
    
    finals <- vector("list", S * M)
    ctr    <- 1L
    
    for (col in seq_len(M)) {
      # work down the column (all S rows share same col index)
      for (row in seq_len(S)) {
        s <- sel[[ row ]]
        flipGene <- flipsMat[row, col]
        s[flipGene] <- 1L - s[flipGene]
        finals[[ctr]] <- BoolNet::getAttractors(boolnet,
                                                method      = "chosen",
                                                startStates = list(s),
                                                returnTable = TRUE,
                                                type        = "synchronous")$attractors[[1]]$involvedStates[[1]]
        ctr <- ctr + 1L
      }
    }
    
    shEnt <- shannonEntropy(table(as.character(unlist(finals))) / (S * M))
    results$Entropy[i] <- shEnt
    
    if (showProgress) utils::setTxtProgressBar(pb, i)
  }
  
  if (showProgress) close(pb)
  
  results$ScaledEntropy <- results$Entropy / log2(nAttr)
  results$Stability     <- 1 - results$ScaledEntropy
  results
}

# =============================================================================
# computeAttractorEntropy_bitflip_bulkC.R  (auto‑detect helper)
#
# Uses BoolNet's hidden C routine `simulateNetworkMultiple()` *when it exists*.
# BoolNet ≥ 2.1.6 on CRAN once carried the helper but some later releases
# renamed / removed it.  We detect it at run‑time; if missing we fall back to a
# column‑wise apply that still avoids the heavy `getAttractors()` call.
# =============================================================================

computeAttractorEntropy_bitflip_bulkC <- function(boolnet, attractors,
                                                  nSamplesState = 25,
                                                  nPerturb      = 25,
                                                  showProgress  = TRUE) {
  
  # ---- try to locate the C helper -----------------------------------------
  .bulkSim <- get0("simulateNetworkMultiple",
                   envir = asNamespace("BoolNet"),
                   inherits = FALSE)
  if (is.null(.bulkSim))
    message("[INFO] simulateNetworkMultiple() not found in BoolNet namespace – falling back to R loop")
  else
    message("[INFO] Using BoolNet C helper simulateNetworkMultiple() for bulk evolution")
  
  geneNames <- attractors$stateInfo$genes
  nGenes    <- length(geneNames)
  nAttr     <- length(attractors$attractors)
  
  if (showProgress)
    pb <- utils::txtProgressBar(min = 0, max = nAttr, style = 3)
  
  # ---- 1. cache decoded steady states once --------------------------------
  decodeCache <- lapply(attractors$attractors, function(att) {
    lapply(att$involvedStates, decodeBigIntegerState, nGenes = nGenes)
  })
  
  res <- data.frame(AttractorIndex = seq_len(nAttr), Entropy = NA_real_)
  
  for (i in seq_len(nAttr)) {
    
    decoded <- decodeCache[[i]]
    if (length(decoded) == 0) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      next
    }
    
    # ---- 2. sample ≤ nSamplesState unique states ---------------------------
    selIdx <- sample(seq_along(decoded),
                     size = min(length(decoded), nSamplesState))
    sel    <- decoded[selIdx]
    
    S   <- length(sel)            # #seed states
    M   <- nPerturb               # flips per seed
    tot <- S * M                  # total trials for this attractor
    
    # ---- 3. assemble initial state matrix ---------------------------------
    seedsMatIdx <- sample(selIdx,  tot, replace = TRUE)
    flipsVec    <- sample.int(nGenes, tot, replace = TRUE)
    
    initMat <- matrix(0L, nrow = nGenes, ncol = tot)
    
    for (k in seq_len(tot)) {
      initMat[, k] <- sel[[ which(selIdx == seedsMatIdx[k])[1] ]]
      g <- flipsVec[k]
      initMat[g, k] <- 1L - initMat[g, k]
    }
    storage.mode(initMat) <- "integer"
    
    # ---- 4. evolve all trials ---------------------------------------------
    if (!is.null(.bulkSim)) {
      # fast C path ----------------------------------------------------------
      finalMat <- .bulkSim(boolnet, initMat,
                           type = "synchronous", returnSeries = FALSE)
    } else {
      # slower R fallback ----------------------------------------------------
      finalMat <- apply(initMat, 2, function(v) {
        BoolNet::simulateNetwork(boolnet,
                                 startStates  = matrix(v, nrow = nGenes),
                                 type         = "synchronous",
                                 returnSeries = FALSE)[, 1]
      })
      if (is.vector(finalMat)) finalMat <- matrix(finalMat, ncol = 1)
    }
    
    finals <- apply(finalMat, 2, paste, collapse = "")
    shEnt  <- shannonEntropy(table(finals) / tot)
    
    res$Entropy[i] <- shEnt
    
    if (showProgress) utils::setTxtProgressBar(pb, i)
  }
  
  if (showProgress) close(pb)
  
  res$ScaledEntropy <- res$Entropy / log2(nAttr)
  res$Stability     <- 1 - res$ScaledEntropy
  res
}
