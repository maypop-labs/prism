# =============================================================================
# 08 - Attractor Analysis.R
# Purpose: Complete Boolean network analysis and rejuvenation target identification
# =============================================================================

# Load required packages and initialize paths
source("managers/attractorManager.R")
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")
config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)

message("Starting comprehensive attractor analysis...")
message("Loading trajectory data and Boolean network...")

# Load all required data
cds <- loadMonocle3(ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
geneMap <- loadObject(ptPaths$geneMap, config, "gene mapping")
attractors <- loadObject(ptPaths$attractors, config, "attractors")
boolnet <- loadObject(ptPaths$boolNet, config, "BoolNet network")

# =============================================================================
# PART 1: Age Score Computation
# =============================================================================

message("Part 1: Computing age scores for ", length(attractors$attractors), " attractors...")

# Extract binary expression and pseudotime
matBin <- assay(cds, "binary")
colData(cds)$Pseudotime <- pseudotime(cds)
pseudotimeValues <- colData(cds)$Pseudotime

# Order cells by pseudotime for young/old reference computation
cellOrder <- order(pseudotimeValues, na.last = NA)
nCells    <- length(cellOrder)
q20       <- floor(0.2 * nCells)

# Compute reference vectors from trajectory extremes
youngVec <- rowMeans(matBin[, cellOrder[1:q20]], na.rm = TRUE)
oldVec   <- rowMeans(matBin[, cellOrder[(nCells - q20 + 1):nCells]], na.rm = TRUE)

# Map attractor genes to expression data
sanitizedNames <- attractors$stateInfo$genes
originalNames <- geneMap$OriginalName[match(sanitizedNames, geneMap$SanitizedName)]
netGenes <- intersect(originalNames[!is.na(originalNames)], rownames(matBin))

# Subset reference vectors to valid genes
youngVec <- youngVec[netGenes]
oldVec   <- oldVec[netGenes]

# Compute attractor age scores
attractorList   <- attractors$attractors
nAttractors     <- length(attractorList)
attractorScores <- numeric(nAttractors)

for (i in seq_len(nAttractors)) {
  attractor <- attractorList[[i]]
  encodedState <- attractor$involvedStates[[1]]
  decoded      <- decodeBigIntegerState(encodedState, length(sanitizedNames))
  names(decoded) <- originalNames
  decoded <- decoded[netGenes]
  
  # Compute age score via projection onto young-old axis
  numerator   <- sum((decoded - youngVec) * (oldVec - youngVec), na.rm = TRUE)
  denominator <- sum((oldVec - youngVec)^2, na.rm = TRUE)
  
  attractorScores[i] <- if (denominator > 1e-10) numerator / denominator else 0
}

# Extract basin sizes
basinSizes <- sapply(attractorList, function(x) {
  if (is.null(x$basinSize)) 1 else x$basinSize
})

# Create results data frame
attractorDf <- data.frame(
  Attractor       = seq_len(nAttractors),
  AgeScore        = attractorScores,
  BasinSize       = basinSizes / sum(basinSizes),
  AttractorSize   = sapply(attractorList, function(x) length(x$involvedStates)),
  stringsAsFactors = FALSE
)

# Create reference vectors object
referenceVectors <- list(
  youngVec = youngVec,
  oldVec = oldVec,
  netGenes = netGenes
)

message("Part 1 completed: Age scores computed for all attractors")

# =============================================================================
# PART 2: Entropy and Stability Analysis
# =============================================================================

message("Part 2: Computing entropy and stability metrics...")

# Configure entropy analysis parameters
nSamplesState <- config$nSamplesState %||% 15
nPerturb      <- config$nPerturb %||% 25

# Compute attractor entropy and stability
entropyDf <- computeAttractorEntropy(
  boolnet       = boolnet,
  attractors    = attractors,
  nSamplesState = nSamplesState,
  nPerturb      = nPerturb
)

# Merge entropy data with age scores
attractorDfScores <- merge(
  x    = entropyDf,
  y    = attractorDf,
  by.x = "AttractorIndex",
  by.y = "Attractor",
  all.x = TRUE
)

# Calculate composite aging scores
attractorDfScores$AttractorScore <- with(attractorDfScores, {
  stability <- ifelse(is.na(Stability), 0.5, Stability)
  basinSize <- ifelse(is.na(BasinSize), 0, BasinSize)
  ageScore  <- ifelse(is.na(AgeScore), 0.5, AgeScore)
  
  stability * basinSize * ageScore
})

# Compute overall network aging score
overallAgingScore <- sum(attractorDfScores$AttractorScore, na.rm = TRUE)
initialScore <- overallAgingScore

# Create network summary
networkSummary <- data.frame(
  OverallAgingScore = overallAgingScore,
  TotalAttractors = nrow(attractorDfScores),
  MeanStability = mean(attractorDfScores$Stability, na.rm = TRUE),
  stringsAsFactors = FALSE
)

message("Part 2 completed: Network aging score = ", round(overallAgingScore, 4))

# =============================================================================
# PART 3: Perturbation Analysis
# =============================================================================

message("Part 3: Starting systematic perturbation analysis...")

# Check if network is suitable for perturbation
minThreshold <- config$attInitialThreshold %||% 0.3
if (initialScore <= minThreshold) {
  stop("Initial aging score (", round(initialScore, 4), 
       ") below threshold (", minThreshold, "). Network not suitable for perturbation.")
}

message("Initial aging score: ", round(initialScore, 4), " - proceeding with perturbation")

# Robust perturbation scoring function
scorePerturbationRobust <- function(network, geneNames, values, youngVec, oldVec) {
  perturbedNet <- fixGenes(network, fixIndices = geneNames, values = values)
  
  attractorResults <- getAttractors(perturbedNet, method = "random", 
                                   startStates = 100, type = "synchronous")
  
  if (length(attractorResults$attractors) == 0) {
    return(list(score = NA_real_, success = FALSE))
  }
  
  # Compute weighted aging score
  ageScores <- numeric(length(attractorResults$attractors))
  basinSizes <- numeric(length(attractorResults$attractors))
  
  for (i in seq_along(attractorResults$attractors)) {
    attractor <- attractorResults$attractors[[i]]
    
    if (length(attractor$involvedStates) > 0) {
      decoded <- decodeBigIntegerState(attractor$involvedStates[[1]], 
                                      length(network$genes))
      names(decoded) <- network$genes
      
      commonGenes <- intersect(names(decoded), names(youngVec))
      if (length(commonGenes) > 0) {
        decodedSub <- decoded[commonGenes]
        youngSub <- youngVec[commonGenes]
        oldSub <- oldVec[commonGenes]
        
        numerator <- sum((decodedSub - youngSub) * (oldSub - youngSub), na.rm = TRUE)
        denominator <- sum((oldSub - youngSub)^2, na.rm = TRUE)
        
        ageScores[i] <- if (denominator > 1e-10) numerator / denominator else 0.5
      } else {
        ageScores[i] <- 0.5
      }
      
      basinSizes[i] <- attractor$basinSize %||% 1
    }
  }
  
  if (sum(basinSizes) > 0) {
    basinFractions <- basinSizes / sum(basinSizes)
    weightedScore <- sum(basinFractions * ageScores, na.rm = TRUE)
  } else {
    weightedScore <- mean(ageScores, na.rm = TRUE)
  }
  
  return(list(score = weightedScore, success = TRUE))
}

# Single gene perturbation analysis
testGenes <- boolnet$genes
nGenes <- length(testGenes)

message("Testing ", nGenes, " genes for knockdown and overexpression effects...")

# Initialize results
singleKD <- data.frame(
  Gene = testGenes,
  AgingScore = NA_real_,
  Delta = NA_real_,
  Success = FALSE,
  stringsAsFactors = FALSE
)

singleOE <- data.frame(
  Gene = testGenes,
  AgingScore = NA_real_,
  Delta = NA_real_,
  Success = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(testGenes)) {
  gene <- testGenes[i]
  
  # Test knockdown
  kdResult <- scorePerturbationRobust(boolnet, gene, 0, youngVec, oldVec)
  singleKD$AgingScore[i] <- kdResult$score
  singleKD$Success[i] <- kdResult$success
  singleKD$Delta[i] <- kdResult$score - initialScore
  
  # Test overexpression
  oeResult <- scorePerturbationRobust(boolnet, gene, 1, youngVec, oldVec)
  singleOE$AgingScore[i] <- oeResult$score
  singleOE$Success[i] <- oeResult$success
  singleOE$Delta[i] <- oeResult$score - initialScore
}

# Results summary
validKD <- sum(singleKD$Success)
validOE <- sum(singleOE$Success)
kdImprovements <- sum(singleKD$Delta < -0.01, na.rm = TRUE)
oeImprovements <- sum(singleOE$Delta < -0.01, na.rm = TRUE)

message("Part 3 completed - Perturbation analysis results:")
message("  - Valid KD results: ", validKD, "/", nGenes)
message("  - Valid OE results: ", validOE, "/", nGenes)
message("  - KD improvements: ", kdImprovements)
message("  - OE improvements: ", oeImprovements)

# Sort by improvement
singleKD <- singleKD[order(singleKD$AgingScore), ]
singleOE <- singleOE[order(singleOE$AgingScore), ]

# Create perturbation analysis summary
analysisSummary <- data.frame(
  InitialAgingScore = initialScore,
  GenesAnalyzed = nGenes,
  ValidKDResults = validKD,
  ValidOEResults = validOE,
  KDImprovements = kdImprovements,
  OEImprovements = oeImprovements,
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 4: Target Ranking and Conclusion
# =============================================================================

message("Part 4: Consolidating results and ranking rejuvenation targets...")

# Combine and analyze results
singleKD$Mode <- "Knockdown"
singleOE$Mode <- "Overexpression"
allResults <- rbind(singleKD, singleOE)

# Filter for successful results and improvements
successful <- allResults[allResults$Success & !is.na(allResults$Delta), ]
improvements <- successful[successful$Delta < -0.01, ]

message("Target identification results:")
message("  - Total tests: ", nrow(allResults))
message("  - Successful: ", nrow(successful))
message("  - Significant improvements: ", nrow(improvements))

# Initialize empty objects for saving
pivotTable <- data.frame()
topTargetsResults <- data.frame()

if (nrow(improvements) > 0) {
  # Sort by improvement (most negative delta first)
  improvements <- improvements[order(improvements$Delta), ]
  
  # Add ranking
  improvements$Rank <- 1:nrow(improvements)
  improvements$PercentImprovement <- abs(improvements$Delta) * 100 / abs(improvements$Delta[1])
  
  message("\nTop 5 rejuvenation targets:")
  topTargetsResults <- head(improvements[, c("Gene", "Mode", "Delta", "Rank")], 5)
  print(topTargetsResults)
  
  # Create simplified pivot table for TSV export
  pivotData <- improvements[, c("Gene", "Mode", "Delta", "AgingScore")]
  
  # Reshape to wide format manually
  kdData <- pivotData[pivotData$Mode == "Knockdown", c("Gene", "Delta", "AgingScore")]
  oeData <- pivotData[pivotData$Mode == "Overexpression", c("Gene", "Delta", "AgingScore")]
  
  colnames(kdData) <- c("Gene", "KD_Delta", "KD_AgingScore")
  colnames(oeData) <- c("Gene", "OE_Delta", "OE_AgingScore")
  
  pivotTable <- merge(kdData, oeData, by = "Gene", all = TRUE)
  
  # Calculate best result per gene
  pivotTable$BestDelta <- pmin(pivotTable$KD_Delta, pivotTable$OE_Delta, na.rm = TRUE)
  pivotTable$BestMode <- ifelse(is.na(pivotTable$KD_Delta), "OE",
                         ifelse(is.na(pivotTable$OE_Delta), "KD",
                         ifelse(pivotTable$KD_Delta < pivotTable$OE_Delta, "KD", "OE")))
  
  # Sort by best result
  pivotTable <- pivotTable[order(pivotTable$BestDelta), ]
  
  message("Part 4 completed: ", nrow(pivotTable), " prioritized targets identified")
  
} else {
  message("Part 4 completed: No significant improvements found")
}

# =============================================================================
# PART 5: Double Gene Perturbation Analysis
# =============================================================================

message("Part 5: Starting double gene perturbation analysis...")

# Configure double perturbation parameters
nTopTargets <- config$nTopTargetsDouble %||% 5
nStatesDouble <- config$nStatesDouble %||% 100

# Get top targets for combination analysis
topKDGenes <- head(singleKD[singleKD$Success & !is.na(singleKD$Delta), ]$Gene, nTopTargets)
topOEGenes <- head(singleOE[singleOE$Success & !is.na(singleOE$Delta), ]$Gene, nTopTargets)

message("Testing combinations of top ", nTopTargets, " targets...")
message("  - Top KD genes: ", paste(topKDGenes, collapse = ", "))
message("  - Top OE genes: ", paste(topOEGenes, collapse = ", "))

# Double perturbation scoring function
scoreDoublePerturbation <- function(network, gene1, gene2, values, youngVec, oldVec) {
  result <- scorePerturbationRobust(network, c(gene1, gene2), values, youngVec, oldVec)
  return(result$score)
}

# Initialize double perturbation results
doubleKD <- data.frame(
  Gene1 = character(),
  Gene2 = character(),
  AgingScore = numeric(),
  Delta = numeric(),
  stringsAsFactors = FALSE
)

doubleOE <- data.frame(
  Gene1 = character(),
  Gene2 = character(),
  AgingScore = numeric(),
  Delta = numeric(),
  stringsAsFactors = FALSE
)

doubleMix <- data.frame(
  KD_Gene = character(),
  OE_Gene = character(),
  AgingScore = numeric(),
  Delta = numeric(),
  stringsAsFactors = FALSE
)

# Double knockdown combinations
if (length(topKDGenes) >= 2) {
  message("Testing double knockdown combinations...")
  
  for (i in 1:(length(topKDGenes) - 1)) {
    for (j in (i + 1):length(topKDGenes)) {
      gene1 <- topKDGenes[i]
      gene2 <- topKDGenes[j]
      
      score <- scoreDoublePerturbation(boolnet, gene1, gene2, c(0, 0), youngVec, oldVec)
      
      doubleKD <- rbind(doubleKD, data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        AgingScore = score,
        Delta = score - initialScore,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Double overexpression combinations
if (length(topOEGenes) >= 2) {
  message("Testing double overexpression combinations...")
  
  for (i in 1:(length(topOEGenes) - 1)) {
    for (j in (i + 1):length(topOEGenes)) {
      gene1 <- topOEGenes[i]
      gene2 <- topOEGenes[j]
      
      score <- scoreDoublePerturbation(boolnet, gene1, gene2, c(1, 1), youngVec, oldVec)
      
      doubleOE <- rbind(doubleOE, data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        AgingScore = score,
        Delta = score - initialScore,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Mixed perturbations (KD + OE)
if (length(topKDGenes) >= 1 && length(topOEGenes) >= 1) {
  message("Testing mixed KD/OE combinations...")
  
  for (kdGene in topKDGenes) {
    for (oeGene in topOEGenes) {
      score <- scoreDoublePerturbation(boolnet, kdGene, oeGene, c(0, 1), youngVec, oldVec)
      
      doubleMix <- rbind(doubleMix, data.frame(
        KD_Gene = kdGene,
        OE_Gene = oeGene,
        AgingScore = score,
        Delta = score - initialScore,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sort results by improvement
doubleKD <- doubleKD[order(doubleKD$AgingScore), ]
doubleOE <- doubleOE[order(doubleOE$AgingScore), ]
doubleMix <- doubleMix[order(doubleMix$AgingScore), ]

# Count significant improvements
doubleKDImprovements <- sum(doubleKD$Delta < -0.01, na.rm = TRUE)
doubleOEImprovements <- sum(doubleOE$Delta < -0.01, na.rm = TRUE)
doubleMixImprovements <- sum(doubleMix$Delta < -0.01, na.rm = TRUE)

message("Part 5 completed - Double perturbation results:")
message("  - Double KD combinations tested: ", nrow(doubleKD))
message("  - Double OE combinations tested: ", nrow(doubleOE))
message("  - Mixed KD/OE combinations tested: ", nrow(doubleMix))
message("  - Double KD improvements: ", doubleKDImprovements)
message("  - Double OE improvements: ", doubleOEImprovements)
message("  - Mixed improvements: ", doubleMixImprovements)

# Find best double perturbation overall
bestDoubleScore <- NA_real_
bestDoubleDelta <- NA_real_
bestDoubleType <- "None"
bestDoubleGenes <- "None"

if (nrow(doubleKD) > 0 && !is.na(doubleKD$AgingScore[1])) {
  bestDoubleScore <- doubleKD$AgingScore[1]
  bestDoubleDelta <- doubleKD$Delta[1]
  bestDoubleType <- "Double KD"
  bestDoubleGenes <- paste(doubleKD$Gene1[1], "+", doubleKD$Gene2[1])
}

if (nrow(doubleOE) > 0 && !is.na(doubleOE$AgingScore[1]) && 
    (is.na(bestDoubleScore) || doubleOE$AgingScore[1] < bestDoubleScore)) {
  bestDoubleScore <- doubleOE$AgingScore[1]
  bestDoubleDelta <- doubleOE$Delta[1]
  bestDoubleType <- "Double OE"
  bestDoubleGenes <- paste(doubleOE$Gene1[1], "+", doubleOE$Gene2[1])
}

if (nrow(doubleMix) > 0 && !is.na(doubleMix$AgingScore[1]) && 
    (is.na(bestDoubleScore) || doubleMix$AgingScore[1] < bestDoubleScore)) {
  bestDoubleScore <- doubleMix$AgingScore[1]
  bestDoubleDelta <- doubleMix$Delta[1]
  bestDoubleType <- "Mixed KD/OE"
  bestDoubleGenes <- paste(doubleMix$KD_Gene[1], "(KD) +", doubleMix$OE_Gene[1], "(OE)")
}

if (!is.na(bestDoubleDelta)) {
  message("\nBest double perturbation:")
  message("  Type: ", bestDoubleType)
  message("  Genes: ", bestDoubleGenes)
  message("  Improvement: ", round(bestDoubleDelta, 4))
}

# =============================================================================
# Save All Results
# =============================================================================

if (config$saveResults) {
  message("Saving all analysis results...")
  
  # Part 1 outputs
  saveObject(attractorDf, ptPaths$attractorDf, config, "attractor age scores")
  saveObject(referenceVectors, ptPaths$referenceVectors, config, "reference vectors")
  
  # Part 2 outputs  
  saveObject(attractorDfScores, ptPaths$attractorScoresCombined, config, "combined attractor scores")
  saveObject(networkSummary, ptPaths$networkAgingSummary, config, "network aging summary")
  
  # Part 3 outputs
  saveObject(singleKD, ptPaths$singleTargetsKD, config, "single gene KD results")
  saveObject(singleOE, ptPaths$singleTargetsOE, config, "single gene OE results")
  saveObject(analysisSummary, ptPaths$perturbationSummary, config, "perturbation summary")
  
  # Part 4 outputs
  if (nrow(improvements) > 0) {
    saveObject(improvements, ptPaths$finalTargetRankings, config, "final target rankings")
    
    # Save pivot table for easy viewing
    if (nrow(pivotTable) > 0) {
      write.table(pivotTable, ptPaths$targetSummaryTsv, sep = "\t", row.names = FALSE, quote = FALSE)
      message("Target summary saved: ", nrow(pivotTable), " prioritized targets")
    }
  } else {
    message("No targets to save")
  }
  
  # Part 5 outputs
  if (nrow(doubleKD) > 0) {
    saveObject(doubleKD, ptPaths$doubleTargetsKD, config, "double KD results")
  }
  if (nrow(doubleOE) > 0) {
    saveObject(doubleOE, ptPaths$doubleTargetsOE, config, "double OE results")
  }
  if (nrow(doubleMix) > 0) {
    saveObject(doubleMix, ptPaths$doubleTargetsMix, config, "mixed KD/OE results")
  }
  
  # Generate comprehensive all-targets summary TSV
  message("Generating comprehensive target summary...")
  
  # Initialize comprehensive results list
  allTargetsResults <- data.frame(
    Rank = integer(),
    Target = character(),
    Type = character(),
    AgingScore = numeric(),
    Delta = numeric(),
    PercentImprovement = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add single gene results
  if (nrow(singleKD) > 0) {
    singleKDFormatted <- data.frame(
      Target = singleKD$Gene,
      Type = "Single KD",
      AgingScore = singleKD$AgingScore,
      Delta = singleKD$Delta,
      stringsAsFactors = FALSE
    )
    singleKDFormatted <- singleKDFormatted[singleKD$Success & !is.na(singleKD$Delta), ]
    allTargetsResults <- rbind(allTargetsResults, singleKDFormatted)
  }
  
  if (nrow(singleOE) > 0) {
    singleOEFormatted <- data.frame(
      Target = singleOE$Gene,
      Type = "Single OE",
      AgingScore = singleOE$AgingScore,
      Delta = singleOE$Delta,
      stringsAsFactors = FALSE
    )
    singleOEFormatted <- singleOEFormatted[singleOE$Success & !is.na(singleOE$Delta), ]
    allTargetsResults <- rbind(allTargetsResults, singleOEFormatted)
  }
  
  # Add double KD results
  if (nrow(doubleKD) > 0) {
    doubleKDFormatted <- data.frame(
      Target = paste(doubleKD$Gene1, "+", doubleKD$Gene2),
      Type = "Double KD",
      AgingScore = doubleKD$AgingScore,
      Delta = doubleKD$Delta,
      stringsAsFactors = FALSE
    )
    doubleKDFormatted <- doubleKDFormatted[!is.na(doubleKD$Delta), ]
    allTargetsResults <- rbind(allTargetsResults, doubleKDFormatted)
  }
  
  # Add double OE results
  if (nrow(doubleOE) > 0) {
    doubleOEFormatted <- data.frame(
      Target = paste(doubleOE$Gene1, "+", doubleOE$Gene2),
      Type = "Double OE",
      AgingScore = doubleOE$AgingScore,
      Delta = doubleOE$Delta,
      stringsAsFactors = FALSE
    )
    doubleOEFormatted <- doubleOEFormatted[!is.na(doubleOE$Delta), ]
    allTargetsResults <- rbind(allTargetsResults, doubleOEFormatted)
  }
  
  # Add mixed KD/OE results
  if (nrow(doubleMix) > 0) {
    doubleMixFormatted <- data.frame(
      Target = paste(doubleMix$KD_Gene, "(KD) +", doubleMix$OE_Gene, "(OE)"),
      Type = "Mixed KD/OE",
      AgingScore = doubleMix$AgingScore,
      Delta = doubleMix$Delta,
      stringsAsFactors = FALSE
    )
    doubleMixFormatted <- doubleMixFormatted[!is.na(doubleMix$Delta), ]
    allTargetsResults <- rbind(allTargetsResults, doubleMixFormatted)
  }
  
  # Sort by improvement (best delta first) and add rankings
  if (nrow(allTargetsResults) > 0) {
    allTargetsResults <- allTargetsResults[order(allTargetsResults$Delta), ]
    allTargetsResults$Rank <- 1:nrow(allTargetsResults)
    
    # Calculate percent improvement relative to best result
    if (allTargetsResults$Delta[1] < 0) {
      allTargetsResults$PercentImprovement <- abs(allTargetsResults$Delta) * 100 / abs(allTargetsResults$Delta[1])
    } else {
      allTargetsResults$PercentImprovement <- 0
    }
    
    # Save comprehensive summary
    allTargetsPath <- gsub("_targetsummary", "_alltargets", ptPaths$targetSummaryTsv)
    write.table(allTargetsResults, allTargetsPath, sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("Comprehensive target summary saved: ", nrow(allTargetsResults), " total targets")
    message("  - File: ", basename(allTargetsPath))
    message("  - Best target: ", allTargetsResults$Target[1], " (", allTargetsResults$Type[1], ")")
    message("  - Best improvement: ", round(allTargetsResults$Delta[1], 4))
    
    # Also save top 20 for quick review
    top20Targets <- head(allTargetsResults, 20)
    top20Path <- gsub("_targetsummary", "_top20targets", ptPaths$targetSummaryTsv)
    write.table(top20Targets, top20Path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  - Top 20 saved: ", basename(top20Path))
    
  } else {
    message("No targets found for comprehensive summary")
  }
}

# =============================================================================
# Final Summary
# =============================================================================

message("\n" , paste(rep("=", 60), collapse = ""))
message("PRISM ATTRACTOR ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("Network aging score: ", round(overallAgingScore, 4))
message("Genes analyzed: ", nGenes)
message("Successful perturbations: ", nrow(successful), "/", nrow(allResults))

if (nrow(improvements) > 0) {
  bestKD <- min(singleKD$Delta, na.rm = TRUE)
  bestOE <- min(singleOE$Delta, na.rm = TRUE)
  bestSingle <- min(bestKD, bestOE)
  
  message("Single gene improvements found: ", nrow(improvements))
  message("Best single improvement: ", round(bestSingle, 4), " (", round(abs(bestSingle) * 100, 2), "% reduction)")
  
  if (nrow(topTargetsResults) > 0) {
    message("\nTop single target: ", topTargetsResults$Gene[1], " (", topTargetsResults$Mode[1], ")")
    message("  Improvement: ", round(topTargetsResults$Delta[1], 4))
  }
} else {
  message("No single gene improvements identified")
  bestSingle <- NA_real_
}

# Double perturbation summary
totalDoubleCombinations <- nrow(doubleKD) + nrow(doubleOE) + nrow(doubleMix)
totalDoubleImprovements <- doubleKDImprovements + doubleOEImprovements + doubleMixImprovements

if (totalDoubleCombinations > 0) {
  message("\nDouble perturbation combinations tested: ", totalDoubleCombinations)
  message("Double perturbation improvements found: ", totalDoubleImprovements)
  
  if (!is.na(bestDoubleDelta)) {
    message("Best double improvement: ", round(bestDoubleDelta, 4), " (", round(abs(bestDoubleDelta) * 100, 2), "% reduction)")
    message("\nTop double target: ", bestDoubleGenes, " (", bestDoubleType, ")")
    
    # Compare single vs double
    if (!is.na(bestSingle) && !is.na(bestDoubleDelta)) {
      if (bestDoubleDelta < bestSingle) {
        synergy <- round(abs(bestDoubleDelta - bestSingle), 4)
        message("\nSynergy detected: ", synergy, " additional improvement from combination")
      } else {
        message("\nNo synergy: Single targets perform better than combinations")
      }
    }
  } else {
    message("No significant double perturbation improvements found")
  }
} else {
  message("\nNo double perturbations tested")
}

message(paste(rep("=", 60), collapse = ""))
