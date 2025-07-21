# =============================================================================
# 05 - Smoothing
#
# Smooth expression profiles along pseudotime using a moving average window.
# Store smoothed values as a new assay and save the updated Monocle3 object.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(monocle3)
library(Seurat)
library(Matrix)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path   <- paste0(config$rootPath, "results/monocle3/")
plotPath       <- paste0(config$rootPath, "results/plots/")
rdsPath        <- paste0(config$rootPath, "results/rds/")
tsvPath        <- paste0(config$rootPath, "results/tsv/")
txtPath        <- paste0(config$rootPath, "results/txt/")
cellTypes      <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType       <- showCellTypeMenu(cellTypes)
trajNamesFile  <- readRDS(paste0(rdsPath, "retained_trajectories_", cellType, ".rds"))
cellTrajectory <- showTrajectoryMenu(trajNamesFile)
cdsPath        <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory)
savePath       <- paste0(cdsPath, "_smoothed")

dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,  recursive = TRUE, showWarnings = FALSE)

# --- Load Monocle3 Object ---
if (!dir.exists(cdsPath)) stop("Monocle3 directory not found: ", cdsPath)
message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)

cat("\014")
cat("\n")

# --- Ensure Pseudotime is Computed ---
if (is.null(colData(cds)$Pseudotime)) {
  message("Computing pseudotime")
  colData(cds)$Pseudotime <- pseudotime(cds)
}

# --- Moving Average Smoothing ---
cellOrder <- order(colData(cds)$Pseudotime)
nCells    <- length(cellOrder)
winSize   <- max(5, floor(config$winSizePercent * nCells))
halfWin   <- floor(winSize / 2)
exprMat   <- assay(cds, "counts")
smoothMat <- matrix(0, nrow = nrow(exprMat), ncol = ncol(exprMat), dimnames = dimnames(exprMat))
message("Smoothing expression values with window size: ", winSize)
for (i in seq_len(nCells)) {
  idx <- cellOrder[max(1, i - halfWin):min(nCells, i + halfWin)]
  smoothMat[, cellOrder[i]] <- rowMeans(exprMat[, idx, drop = FALSE])
}
assay(cds, "smoothed_expr") <- smoothMat

cat("\014")
cat("\n")

# --- Save Results ---
if (config$saveResults) {
  message("Saving smoothed Monocle3 object to: ", savePath)
  save_monocle_objects(cds = cds, directory_path = savePath, comment = cellType)
}

message("Done!")
