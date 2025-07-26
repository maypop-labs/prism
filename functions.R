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
