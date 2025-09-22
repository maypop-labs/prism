# =============================================================================
# 03 - Pseudotime
#
# Convert a cell-type-specific Seurat object to a Monocle3 CellDataSet,
# perform pseudotime trajectory inference, identify valid branches,
# and save valid sub-trajectories.
# =============================================================================

# --- Initialization ---
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

# --- Load Data ---
seuratObj <- loadObject(ctPaths$seuratObject, config, "cell type Seurat object")

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
      saveMonocle3(subCds, ptPaths$monocle3, config, "pseudotime trajectory")
    }
  }
}

# ---- Post‑processing of Branch Statistics ---
branch_stats <- branch_stats |>
  dplyr::filter(!is.na(corWithAge) & corWithAge >= 0) |>
  dplyr::arrange(dplyr::desc(corWithAge))

# --- Save Results ---
if (config$saveResults) {
  saveObject(branch_stats, ctPaths$trajectoryCorrelations, config, "branch statistics report")
  saveObject(retainedLeaves, ctPaths$retainedTrajectories, config, "retained trajectory names")
}

# --- Create and Save Plots ---
if (config$saveResults && length(retainedLeaves) > 0) {
  if (config$verbose) { message("Creating age vs pseudotime plots...") }

  # Create the plots
  trajectoryPlots <- createTrajectoryPlots(
    cds            = cds,
    retainedLeaves = retainedLeaves,
    branchStats    = branch_stats,
    basePaths      = paths$base,
    cellType       = cellType,
    config         = config,
    saveIndividual = TRUE
  )
}

if (config$verbose) { message(retained, " valid trajectory branch(es).") }
message("Done!")
