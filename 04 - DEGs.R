# =============================================================================
# 04 - DEGs
#
# Load a Monocle3 trajectory object and identify genes whose expression
# varies significantly along the pseudotime principal graph.
# Save differentially expressed genes (DEGs) for downstream analysis.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(monocle3)
library(Seurat)
library(Matrix)

# --- Source functions ---
source("functions.R")

# --- Config and options ---
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

dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,  recursive = TRUE, showWarnings = FALSE)

cdsPath      <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory)
degSaveFile  <- paste0(rdsPath, cellType, "_", cellTrajectory, "_degs.rds")
degTsvFile   <- paste0(tsvPath,  cellType, "_", cellTrajectory, "_degs.tsv")

# --- Load Monocle3 Object ---
if (!dir.exists(cdsPath)) stop("Monocle3 directory not found: ", cdsPath)
message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)

cat("\014")
cat("\n")

# --- DEG Testing Along Pseudotime ---
message("Performing graph_test for differential gene expression")
degTable <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
message("DEG test complete")

# --- Filter Significant Genes ---
degTable <- degTable[degTable$q_value < config$fdrLevel, ]
degTable <- degTable[order(degTable$q_value), ]

cat("\014")
cat("\n")

# --- Save Output ---
if (config$saveResults) {
  message("Saving DEGs to: ", degSaveFile)
  saveRDS(degTable, file = degSaveFile)
  
  message("Saving DEG report to: ", degTsvFile)
  reportCols <- c("gene_id", "gene_short_name", "status", "p_value", "q_value",
                  "morans_I", "mean_expression", "num_cells_expressed")
  reportCols <- intersect(reportCols, colnames(degTable))
  degReport <- as.data.frame(degTable[, reportCols, drop = FALSE])
  write.table(degReport, file = degTsvFile, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Number of significant DEGs: ", nrow(degTable))
message("Done!")
