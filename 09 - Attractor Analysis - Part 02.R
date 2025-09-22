# =============================================================================
# 09 - Attractor Analysis - Part 02.R
# Purpose: Compute entropy and stability metrics for Boolean network attractors
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

message("Loading Boolean network and attractor age scores...")

# Load data
boolnet <- loadObject(ptPaths$boolNet, config, "BoolNet network")
attractors <- loadObject(ptPaths$attractors, config, "attractors")
attractorDf <- loadObject(ptPaths$attractorDf, config, "attractor age scores")

# Configure entropy analysis parameters
nSamplesState <- config$nSamplesState %||% 15
nPerturb      <- config$nPerturb %||% 25

message("Computing entropy for ", length(attractors$attractors), " attractors...")

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

message("Network aging score: ", round(overallAgingScore, 4))

# Save results
if (config$saveResults) {
  saveObject(attractorDfScores, ptPaths$attractorScoresCombined, config, "combined attractor scores")
  
  networkSummary <- data.frame(
    OverallAgingScore = overallAgingScore,
    TotalAttractors = nrow(attractorDfScores),
    MeanStability = mean(attractorDfScores$Stability, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  saveObject(networkSummary, ptPaths$networkAgingSummary, config, "network aging summary")
}
message("Done!")
