# =============================================================================
# 03 - Pseudotime
#
# Convert a cell-type-specific Seurat object to a Monocle3 CellDataSet,
# perform pseudotime trajectory inference, identify valid branches,
# and save valid sub-trajectories.
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
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
if (config$verbose) { message("Loading Seurat object for cell type: ", cellType) }
seuratObj <- readRDS(ctPaths$seuratObject)

# --- Convert to Monocle3 CellDataSet ---
cds <- convertSeuratToCDS(seuratObj)
rootAge = min(config$ages, na.rm = TRUE)

# --- Pseudotime Rooting ---
if (config$verbose) { message("Running pseudotime analysis from age ", rootAge) }
rootCells <- getRootCells(cds, rootAge)
if (length(rootCells) == 0) { stop("No cells found for root age: ", rootAge) }
cds <- runPseudotime(cds, rootCells, verbose = config$verbose)

# --- Identify Valid Branches ---
if (config$verbose) { message("Searching for trajectory branches with high age correlation") }
trajGraph     <- principal_graph(cds)[["UMAP"]]
closestVertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
leafNodes     <- names(which(degree(trajGraph) == 1))
rootVertices  <- unique(V(trajGraph)$name[closestVertex[rootCells, 1]])
retained      <- 0

# --- create a table for storing leaf node names and correlations ---
branch_stats <- tibble::tibble(
  leaf        = character(),
  corWithAge  = numeric(),
  method      = character(),
  nCells      = numeric(),
  nDonors     = numeric(),
  pValue      = numeric()
)
retainedLeaves <- character()

for (leaf in leafNodes) {
  ptPaths            <- getTrajectoryFilePaths(paths$base, cellType, leaf)
  sp                 <- shortest_paths(trajGraph, from = rootVertices, to = leaf, weights = NA)$vpath[[1]]
  pathNodeIdx        <- as.integer(sub("Y_", "", as_ids(sp)))
  cellBarcodes       <- rownames(closestVertex)[closestVertex[, 1] %in% pathNodeIdx]
  subCds             <- cds[, cellBarcodes]
  pt                 <- pseudotime(subCds)
  nCells             <- length(pt)
  correlation_result <- calculateTrajectoryCorrelation(subCds)
  corWithAge         <- correlation_result$correlation
  pValue             <- correlation_result$p_value
  
  # Record the statistics 
  branch_stats <- dplyr::add_row(branch_stats, 
                                 leaf       = leaf, 
                                 corWithAge = corWithAge,
                                 method     = correlation_result$method,
                                 nCells     = correlation_result$n_cells,
                                 pValue     = pValue,
                                 nDonors    = correlation_result$n_donors
  )
  
  # Apply threshold (using p-value now too)
  if (is.finite(corWithAge) && corWithAge >= config$pseudotimeMinAgeCorrelation && pValue < 0.05) {
    retained <- retained + 1
    retainedLeaves <- c(retainedLeaves, leaf)

    if (config$saveResults) {
      if (config$verbose) { message("Saving trajectory to: ", ptPaths$monocle3) }
      save_monocle_objects(cds = subCds, directory_path = ptPaths$monocle3, comment = cellType)
    }
  }
}

# ---- Post‑processing of branch_stats ---
branch_stats <- branch_stats |>
  dplyr::filter(!is.na(corWithAge) & corWithAge >= 0) |>
  dplyr::arrange(dplyr::desc(corWithAge))

# --- save the trajectory report and the trajectory names vector ---
if (config$saveResults) {
  readr::write_tsv(branch_stats, ctPaths$trajectoryCorrelations)
  if (config$verbose) { message("Saved branch statistics report to ", ctPaths$trajectoryCorrelations) }

  saveRDS(retainedLeaves, ctPaths$retainedTrajectories)
  if (config$verbose) { message("Saved retained trajectory names to ", ctPaths$retainedTrajectories) }
  
}

if (config$saveResults && length(retainedLeaves) > 0) {
  if (config$verbose) { message("Creating age vs pseudotime plots...") }
  
  plotPath <- paths$base$plots
  
  # Create the plots
  trajectoryPlots <- createTrajectoryPlots(
    cds = cds,
    retainedLeaves = retainedLeaves,
    branchStats = branch_stats,
    basePaths = paths$base,
    cellType = cellType,
    config = config,
    saveIndividual = TRUE
  )
  
  if (config$verbose) { message("Plots saved to: ", plotPath) }
}

if (config$verbose) { message(retained, " valid trajectory branch(es).") }
message("Done!")
