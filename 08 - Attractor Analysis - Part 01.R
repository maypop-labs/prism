# =============================================================================
# 09 - Attractor Analysis - Part 01
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

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$seuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path    <- paste0(config$rootPath, "results/monocle3/")
graphMlPath     <- paste0(config$rootPath, "results/graphml/")
plotPath        <- paste0(config$rootPath, "results/plots/")
rdsPath         <- paste0(config$rootPath, "results/rds/")
tsvPath         <- paste0(config$rootPath, "results/tsv/")
txtPath         <- paste0(config$rootPath, "results/txt/")
cellTypes       <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType        <- showCellTypeMenu(cellTypes)
trajNamesFile   <- readRDS(paste0(rdsPath, "retained_trajectories_", cellType, ".rds"))
cellTrajectory  <- showTrajectoryMenu(trajNamesFile)
cdsPath         <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile         <- paste0(rdsPath, cellType, "_", cellTrajectory, "_switch_degs.rds")
edgesFile       <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile       <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile       <- paste0(rdsPath, cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
boolnetFile     <- paste0(rdsPath, cellType, "_", cellTrajectory, "_boolnet.rds")
attractorsFile  <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractors.rds")
attractorDfFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractor_df.rds")
geneMapFile     <- paste0(rdsPath, cellType, "_", cellTrajectory, "_gene_map.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(geneMapFile)) stop("Gene Map RDS file not found: ", geneMapFile)
if (!file.exists(attractorsFile)) stop("Attractor RDS file not found: ", attractorsFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading gene map from: ", geneMapFile)
geneMap <- readRDS(geneMapFile)
message("Loading attractors from: ", attractorsFile)
attractors <- readRDS(attractorsFile)

cat("\014")
cat("\n")

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

cat("\014")
cat("\n")

# --- (Optional) Save ---
if (config$saveResults) {
  message("Saving attractor data frame to: ", attractorDfFile)
  saveRDS(attractorDf, file = attractorDfFile)
}

message("Done!")
