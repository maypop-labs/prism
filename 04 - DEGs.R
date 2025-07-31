# =============================================================================
# 04 - DEGs
#
# Load a Monocle3 trajectory object and identify genes whose expression
# varies significantly along the pseudotime principal graph.
# Save differentially expressed genes (DEGs) for downstream analysis.
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

# --- DEG Testing Along Pseudotime ---
message("Performing graph_test for differential gene expression")
degTable <- graph_test(cds, neighbor_graph = "principal_graph", cores = config$cores)
message("DEG test complete")

# --- Filter Significant Genes ---
degTable <- degTable[degTable$q_value < config$fdrLevel, ]
degTable <- degTable[order(degTable$q_value), ]

cat("\014")
cat("\n")

# --- Save Output ---
if (config$saveResults) {
  message("Saving DEGs to: ", ptPaths$degs)
  saveRDS(degTable, file = ptPaths$degs)
  
  message("Saving DEG report to: ", ptPaths$degsTsv)
  reportCols <- c("gene_id", "gene_short_name", "status", "p_value", "q_value",
                  "morans_I", "mean_expression", "num_cells_expressed")
  reportCols <- intersect(reportCols, colnames(degTable))
  degReport <- as.data.frame(degTable[, reportCols, drop = FALSE])
  write.table(degReport, file = ptPaths$degsTsv, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Number of significant DEGs: ", nrow(degTable))
message("Done!")
