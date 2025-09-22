# =============================================================================
# 11 - Conclusion.R
# Purpose: Consolidate perturbation results and create final target rankings
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

message("Loading perturbation analysis results...")

# Load single gene results
singleKD <- loadObject(ptPaths$singleTargetsKD, config, "single gene KD results")
singleOE <- loadObject(ptPaths$singleTargetsOE, config, "single gene OE results")

# Combine and analyze results
singleKD$Mode <- "Knockdown"
singleOE$Mode <- "Overexpression"
allResults <- rbind(singleKD, singleOE)

# Filter for successful results and improvements
successful <- allResults[allResults$Success & !is.na(allResults$Delta), ]
improvements <- successful[successful$Delta < -0.01, ]

message("Results summary:")
message("  - Total tests: ", nrow(allResults))
message("  - Successful: ", nrow(successful))
message("  - Significant improvements: ", nrow(improvements))

if (nrow(improvements) > 0) {
  # Sort by improvement (most negative delta first)
  improvements <- improvements[order(improvements$Delta), ]
  
  # Add ranking
  improvements$Rank <- 1:nrow(improvements)
  improvements$PercentImprovement <- abs(improvements$Delta) * 100 / abs(improvements$Delta[1])
  
  message("\nTop 5 targets:")
  topTargets <- head(improvements[, c("Gene", "Mode", "Delta", "Rank")], 5)
  print(topTargets)
  
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
  
} else {
  message("No significant improvements found")
  pivotTable <- data.frame()
}

# Save results
if (config$saveResults) {
  # Save detailed results
  saveObject(improvements, ptPaths$finalTargetRankings, config, "final target rankings")
  
  # Save pivot table for easy viewing
  if (nrow(pivotTable) > 0) {
    write.table(pivotTable, ptPaths$targetSummaryTsv, sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("Target analysis complete - saved ", nrow(pivotTable), " prioritized targets")
  } else {
    message("No targets to save")
  }
}

message("Done!")
