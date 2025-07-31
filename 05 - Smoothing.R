# =============================================================================
# 05 - Smoothing
#
# Smooth expression profiles along pseudotime using a moving average window.
# Store smoothed values as a new assay and save the updated Monocle3 object.
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
ensureProjectDirectories(paths)
clearConsole()

# --- Load Monocle3 Object ---
if (!dir.exists(ptPaths$monocle3)) stop("Monocle3 directory not found: ", ptPaths$monocle3)
message("Loading Monocle3 object from: ", ptPaths$monocle3)
cds <- load_monocle_objects(directory_path = ptPaths$monocle3)

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

# --- Save Results ---
if (config$saveResults) {
  message("Saving smoothed Monocle3 object to: ", ptPaths$monocle3Smoothed)
  save_monocle_objects(cds = cds, directory_path = ptPaths$monocle3Smoothed, comment = cellType)
}

message("Done!")
