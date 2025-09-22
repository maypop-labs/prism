# =============================================================================
# 08 - Attractor Analysis - Part 01.R
# Purpose: Compute age scores for Boolean network attractors
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

message("Loading trajectory data and Boolean network attractors...")

# Load data
cds <- loadMonocle3(ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
geneMap <- loadObject(ptPaths$geneMap, config, "gene mapping")
attractors <- loadObject(ptPaths$attractors, config, "attractors")

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

message("Computing age scores for ", length(attractors$attractors), " attractors...")

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

# Save results
if (config$saveResults) {
  saveObject(attractorDf, ptPaths$attractorDf, config, "attractor age scores")
  
  referenceVectors <- list(
    youngVec = youngVec,
    oldVec = oldVec,
    netGenes = netGenes
  )
  
  saveObject(referenceVectors, ptPaths$referenceVectors, config, "reference vectors")
}
message("Done!")
