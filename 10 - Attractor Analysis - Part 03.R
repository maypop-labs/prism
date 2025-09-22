# =============================================================================
# 10 - Attractor Analysis - Part 03.R
# Purpose: Systematic in silico perturbation analysis for rejuvenation targets
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

message("Loading network data for perturbation analysis...")

# Load data
boolnet <- loadObject(ptPaths$boolNet, config, "BoolNet network")

attractorScores <- loadObject(ptPaths$attractorScoresCombined, config, "combined attractor scores")

# Load reference vectors
referenceVectors <- loadObject(ptPaths$referenceVectors, config, "reference vectors")

# Extract data
youngVec <- referenceVectors$youngVec
oldVec <- referenceVectors$oldVec
initialScore <- sum(attractorScores$AttractorScore, na.rm = TRUE)

# Check if network is suitable for perturbation
minThreshold <- config$attInitialThreshold %||% 0.3
if (initialScore <= minThreshold) {
  stop("Initial aging score (", round(initialScore, 4), 
       ") below threshold (", minThreshold, "). Network not suitable for perturbation.")
}

message("Starting perturbation analysis with initial aging score: ", round(initialScore, 4))

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
  
  if (i %% 10 == 0) {
    message("Progress: ", i, "/", nGenes, " genes")
  }
  
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

message("Perturbation analysis completed:")
message("  - Valid KD results: ", validKD, "/", nGenes)
message("  - Valid OE results: ", validOE, "/", nGenes)
message("  - KD improvements: ", kdImprovements)
message("  - OE improvements: ", oeImprovements)

# Sort by improvement
singleKD <- singleKD[order(singleKD$AgingScore), ]
singleOE <- singleOE[order(singleOE$AgingScore), ]

# Save results
if (config$saveResults) {
  saveObject(singleKD, ptPaths$singleTargetsKD, config, "single gene KD results")
  saveObject(singleOE, ptPaths$singleTargetsOE, config, "single gene OE results")
  
  analysisSummary <- data.frame(
    InitialAgingScore = initialScore,
    GenesAnalyzed = nGenes,
    ValidKDResults = validKD,
    ValidOEResults = validOE,
    KDImprovements = kdImprovements,
    OEImprovements = oeImprovements,
    stringsAsFactors = FALSE
  )
  
  saveObject(analysisSummary, ptPaths$perturbationSummary, config, "perturbation summary")

}

message("Done!")
