# =============================================================================
# 11 - Attractor Analysis - Part 01 (Refactored)
#
# Compute an age-associated score for each Boolean attractor by projecting it
# onto a pseudotime-defined trajectory axis from young to old cell states.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(dplyr)
library(gmp)

# --- Source Helpers ---
source("functions.R")

# --- Options ---
options(warn = -1)
options(Seurat.object.assay.version = "v5")
registerDoParallel(cores = 1)
saveResults <- TRUE

# --- Parameters ---
cellType       <- "Keratinocytes"
cellTrajectory <- "Y_447"

cdsPath        <- paste0("E:/datasets/omics/skin/results/monocle3/monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile        <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_switch_degs.rds")
edgesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
geneMapFile    <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_gene_map.rds")
attractorsFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractors.rds")

outputFile     <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractor_df.rds")

# --- Load Data ---
cds         <- load_monocle_objects(directory_path = cdsPath)
degTable    <- readRDS(degFile)
scenicEdges <- readRDS(edgesFile)
gMerged     <- readRDS(graphFile)
boolRules   <- readRDS(rulesFile)
geneMap     <- readRDS(geneMapFile)
attractors  <- readRDS(attractorsFile)

# --- Prepare Binary Expression and Pseudotime ---
matBin <- assay(cds, "binary")
colData(cds)$Pseudotime <- pseudotime(cds)

# --- Compute Young and Old Attractors from Expression ---
cellOrder <- order(colData(cds)$Pseudotime)
nCells    <- length(cellOrder)
q20       <- floor(0.2 * nCells)
youngVec  <- rowMeans(matBin[, cellOrder[1:q20]])
oldVec    <- rowMeans(matBin[, cellOrder[(nCells - q20 + 1):nCells]])

# --- Prepare Gene Set for Comparison ---
sanitizedNames <- attractors$stateInfo$genes
originalNames  <- geneMap$OriginalName[match(sanitizedNames, geneMap$SanitizedName)]

if (any(is.na(originalNames))) stop("NA found in gene mapping.")

netGenes <- intersect(originalNames, rownames(matBin))
youngVec <- youngVec[netGenes]
oldVec   <- oldVec[netGenes]

# --- Compute Attractor Scores ---
attractorList <- attractors$attractors
nAttr         <- length(attractorList)
attractorScores <- numeric(nAttr)

for (i in seq_len(nAttr)) {
  encodedState <- attractorList[[i]]$involvedStates[[1]]
  decoded      <- decodeBigIntegerState(encodedState, length(sanitizedNames))
  names(decoded) <- originalNames
  decoded <- decoded[netGenes]
  
  numerator   <- sum((decoded - youngVec) * (oldVec - youngVec))
  denominator <- sum((oldVec - youngVec)^2)
  score <- numerator / denominator
  
  attractorScores[i] <- score
}

# --- Assemble Output Data Frame ---
basinSizes <- sapply(attractorList, function(x) x$basinSize)
relBasinSizes <- basinSizes / sum(basinSizes)

attractorDf <- data.frame(
  Attractor  = seq_len(nAttr),
  AgeScore   = attractorScores,
  BasinSize  = relBasinSizes
)

# --- (Optional) Save ---
if (saveResults) {
  saveRDS(attractorDf, file = outputFile)
  message("Saved attractor_df to: ", outputFile)
}

message("Done!")
