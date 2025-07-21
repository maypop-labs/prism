# =============================================================================
# 03 - Pseudotime Core
#
# Convert a cell-type-specific Seurat object to a Monocle3 CellDataSet,
# perform pseudotime trajectory inference, identify valid branches,
# and save valid sub-trajectories.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(Matrix)
library(igraph)
library(tools)

# --- Source functions ---
source("functions.R")

# --- Config and options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
trajRootAge <- min(config$ageVec)
plotPath    <- paste0(config$rootPath, "results/plots/")
rdsPath     <- paste0(config$rootPath, "results/rds/")
tsvPath     <- paste0(config$rootPath, "results/tsv/")
txtPath     <- paste0(config$rootPath, "results/txt/")
cellTypes   <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType    <- showCellTypeMenu(cellTypes)
seuratFile  <- paste0(rdsPath, cellType, "_seurat.rds")

dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,  recursive = TRUE, showWarnings = FALSE)

# --- Load Seurat Object ---

if (!file.exists(seuratFile)) stop("Seurat RDS file not found: ", seuratFile)
message("Loading Seurat object for cell type: ", cellType)
seuratObj <- readRDS(seuratFile)

cat("\014")
cat("\n")

# --- Convert to Monocle3 CellDataSet ---
convertSeuratToCDS <- function(seuratObj) {
  rawCounts <- LayerData(seuratObj, layer = "counts")
  genes <- rownames(rawCounts)
  cellMeta <- seuratObj@meta.data
  geneMeta <- data.frame(gene_id = genes, gene_short_name = genes, row.names = genes)
  new_cell_data_set(rawCounts, cellMeta, geneMeta)
}

cds <- convertSeuratToCDS(seuratObj)

# --- Pseudotime Rooting ---
getRootCells <- function(cds, rootAge) {
  meta <- as.data.frame(colData(cds))
  rownames(meta)[meta$age == rootAge]
}

runPseudotime <- function(cds, rootCells) {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "log", scaling = TRUE, verbose = config$verbose)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose = config$verbose)
  cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden", verbose = config$verbose)
  cds <- learn_graph(cds)
  order_cells(cds, reduction_method = "UMAP", root_cells = rootCells)
}

message("Running pseudotime analysis from age ", trajRootAge)
rootCells <- getRootCells(cds, trajRootAge)
cds <- runPseudotime(cds, rootCells)

# --- Identify Valid Branches ---
message("Searching for trajectory branches with high age correlation")
trajGraph <- principal_graph(cds)[["UMAP"]]
closestVertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
leafNodes <- names(which(degree(trajGraph) == 1))
rootVertices <- unique(V(trajGraph)$name[closestVertex[rootCells, 1]])
retained <- 0

# --- create a table for storing leaf node names and correlations ---
branch_stats <- tibble::tibble(
  leaf        = character(),   # branch (vertex) ID
  corWithAge  = numeric()      # Pearson r between pseudotime and age
)
retainedLeaves <- character()

for (leaf in leafNodes) {
  sp           <- shortest_paths(trajGraph, from = rootVertices, to = leaf, weights = NA)$vpath[[1]]
  pathNodeIdx  <- as.integer(sub("Y_", "", as_ids(sp)))
  cellBarcodes <- rownames(closestVertex)[closestVertex[, 1] %in% pathNodeIdx]
  subCds       <- cds[, cellBarcodes]
  pt           <- pseudotime(subCds)
  corWithAge   <- suppressWarnings(cor(pt, colData(subCds)$age))
  
  # Record the statistics regardless of whether the branch is retained
  branch_stats <- dplyr::add_row(branch_stats, leaf = leaf, corWithAge = corWithAge)

  if (is.finite(corWithAge) && corWithAge >= config$ageCorrelation) {
    retained <- retained + 1
    retainedLeaves <- c(retainedLeaves, leaf)
    message("--- Branch retained from root to leaf ", leaf)
    message("    → Correlation with age: ", round(corWithAge, 3))
    message("    → Cells: ", ncol(subCds))

    subCds <- runPseudotime(subCds, rootCells = getRootCells(subCds, trajRootAge))

    if (config$saveResults) {
      saveDir <- paste0(config$rootPath, "results/monocle3/monocle3_", cellType, "_", leaf)
      message("Saving trajectory to: ", saveDir)
      save_monocle_objects(cds = subCds, directory_path = saveDir, comment = cellType)
    }
  }
}

# ---- Post‑processing of branch_stats ---
branch_stats <- branch_stats |>
  dplyr::filter(!is.na(corWithAge) & corWithAge >= 0) |>  # drop NA and negatives
  dplyr::arrange(dplyr::desc(corWithAge))                 # sort highest to lowest

cat("\014")
cat("\n")

# --- save the trajectory report and the trajectory names vector ---
if (config$saveResults) {
  reportFile <- paste0(tsvPath, "trajectory_age_correlations_", cellType, ".tsv")
  readr::write_tsv(branch_stats, reportFile)
  message("Saved branch statistics report to ", reportFile)

  namesFile <- paste0(rdsPath, "retained_trajectories_", cellType, ".rds")
  saveRDS(retainedLeaves, namesFile)
  message("Saved retained trajectory names to ", namesFile)
  
}

message(retained, " valid trajectory branch(es).")
message("Done!")
