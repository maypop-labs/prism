# =============================================================================
# 03 - Pseudotime Core
#
# Convert a cell-type-specific Seurat object to a Monocle3 CellDataSet,
# perform pseudotime trajectory inference, identify valid branches,
# and save valid sub-trajectories.
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

config   <- initializeScript()
pathInfo <- initializeInteractivePaths(needsCellType = TRUE)
paths    <- pathInfo$paths
cellType <- pathInfo$cellType
ctPaths  <- getCellTypeFilePaths(paths$base, cellType)
ensureProjectDirectories(paths)
clearConsole()

# --- Load Seurat Object ---
if (!file.exists(ctPaths$seuratObject)) stop("Seurat RDS file not found: ", ctPaths$seuratObject)
message("Loading Seurat object for cell type: ", cellType)
seuratObj <- readRDS(ctPaths$seuratObject)

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
rootAge = min(config$ageVec, na.rm = TRUE)

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

message("Running pseudotime analysis from age ", rootAge)
rootCells <- getRootCells(cds, rootAge)
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

    subCds <- runPseudotime(subCds, rootCells = getRootCells(subCds, rootAge))

    if (config$saveResults) {
      saveDir <- paste0(paths$base$monocle3, "monocle3_", cellType, "_", leaf)
      message("Saving trajectory to: ", saveDir)
      save_monocle_objects(cds = subCds, directory_path = saveDir, comment = cellType)
    }
  }
}

# ---- Post‑processing of branch_stats ---
branch_stats <- branch_stats |>
  dplyr::filter(!is.na(corWithAge) & corWithAge >= 0) |>  # drop NA and negatives
  dplyr::arrange(dplyr::desc(corWithAge))                 # sort highest to lowest

# --- save the trajectory report and the trajectory names vector ---
if (config$saveResults) {
  readr::write_tsv(branch_stats, ctPaths$trajectoryCorrelations)
  message("Saved branch statistics report to ", ctPaths$trajectoryCorrelations)

  saveRDS(retainedLeaves, ctPaths$retainedTrajectories)
  message("Saved retained trajectory names to ", ctPaths$retainedTrajectories)
  
}

message(retained, " valid trajectory branch(es).")
message("Done!")
