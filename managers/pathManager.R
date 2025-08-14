# =============================================================================
# pathManager.R
# Purpose: Centralized path and filename management for PRISM project
# Dependencies: config.yaml
# =============================================================================

# =============================================================================
# Path Builder Functions
# =============================================================================

#' Build base directory paths
#' @param config Configuration object from yaml::read_yaml("config.yaml")
#' @return Named list of base directory paths
buildBasePaths <- function(config) {
  rootPath <- config$rootPath
  
  list(
    root = rootPath,
    cellranger = paste0(rootPath, "cellranger_counts/"),
    results    = paste0(rootPath, "results/"),
    json       = paste0(rootPath, "results/json/"),
    plots      = paste0(rootPath, "results/plots/"),
    rds        = paste0(rootPath, "results/rds/"),
    tsv        = paste0(rootPath, "results/tsv/"),
    txt        = paste0(rootPath, "results/txt/"),
    monocle3   = paste0(rootPath, "results/monocle3/"),
    graphml    = paste0(rootPath, "results/graphml/")
  )
}

#' Build CellRanger input paths
#' @param config Configuration object
#' @return Vector of CellRanger paths for each sample
buildCellrangerPaths <- function(config) {
  basePath <- paste0(config$rootPath, "cellranger_counts/")
  paste0(basePath, config$donorIDs, "_output/outs/filtered_feature_bc_matrix/")
}

# =============================================================================
# Static File Paths (No Cell Type/Trajectory Dependencies)
# =============================================================================

#' Get static file paths that don't depend on cell type or trajectory
#' @param basePaths List from buildBasePaths()
#' @return Named list of static file paths
getStaticFilePaths <- function(basePaths) {
  list(
    # Main Seurat objects
    seuratMerged = paste0(basePaths$rds, "seurat_merged.rds"),
    cellTypes    = paste0(basePaths$rds, "cell_types.rds"),
    cellTypeFreq = paste0(basePaths$tsv, "cell_types.tsv"),
    
    # Main plots
    pcaAllCellsByAge   = paste0(basePaths$plots, "PCA_all_cells_by_age.png"),
    umapAllCellsByAge  = paste0(basePaths$plots, "UMAP_all_cells_by_age.png"),
    pcaAllCellsByType  = paste0(basePaths$plots, "PCA_all_cells_by_cell_type.png"),
    umapAllCellsByType = paste0(basePaths$plots, "UMAP_all_cells_by_cell_type.png")
  )
}

# =============================================================================
# Cell Type-Specific Paths (No Trajectory Dependencies)
# =============================================================================

#' Build cell type-specific file paths
#' @param basePaths List from buildBasePaths()
#' @param cellType Cell type string
#' @return Named list of cell type-specific paths
getCellTypeFilePaths <- function(basePaths, cellType) {

  list(
    # Seurat objects
    seuratObject = paste0(basePaths$rds, "seurat_", cellType, ".rds"),
    retainedTrajectories = paste0(basePaths$rds, "retained_trajectories_", cellType, ".rds"),
    
    # Plots
    pcaPlot                    = paste0(basePaths$plots, "PCA_", cellType, ".png"),
    umapPlot                   = paste0(basePaths$plots, "UMAP_", cellType, ".png"),
    trajectoryCorrelationsPlot = paste0(basePaths$plots, "trajectory_correlations_", cellType, "_summary.png"),
    
    # Analysis outputs
    trajectoryCorrelations = paste0(basePaths$tsv, "trajectory_age_correlations_", cellType, ".tsv")

  )
}

# =============================================================================
# Cell Type + Trajectory-Specific Paths
# =============================================================================

#' Build paths that depend on both cell type and trajectory
#' @param basePaths List from buildBasePaths()
#' @param cellType Cell type string
#' @param trajectory Trajectory string
#' @return Named list of cell type + trajectory-specific paths
getTrajectoryFilePaths <- function(basePaths, cellType, trajectory) {
  basePrefix <- paste0(cellType, "_", trajectory)
  
  list(
    # Monocle3 directories
    monocle3             = paste0(basePaths$monocle3, "monocle3_", basePrefix, "/"),
    monocle3GeneSwitches = paste0(basePaths$monocle3, "monocle3_", basePrefix, "_geneSwitches/"),
    monocle3Smoothed     = paste0(basePaths$monocle3, "monocle3_", basePrefix, "_smoothed/"),

    # Pseudotime Analysis
    ageVsPseudotimePlot         = paste0(basePaths$plots, basePrefix, "_age_vs_pseudotime.png"),
    ageVsPseudotimeCombinedPlot = paste0(basePaths$plots, basePrefix, "_age_vs_pseudotime_combined.png"),
    
    # Analysis intermediate files
    geneSwitches    = paste0(basePaths$rds, basePrefix, "_geneSwitches.rds"),
    geneSwitchesTsv = paste0(basePaths$tsv, basePrefix, "_geneSwitches.tsv"),
    grn             = paste0(basePaths$rds, basePrefix, "_GRN.rds"),
    grnEdges        = paste0(basePaths$rds, basePrefix, "_GRN_edges.rds"),
    grnPreprocessed = paste0(basePaths$rds, basePrefix, "_GRN_preprocessed.rds"),
    
    # Boolean network files
    booleanRules = paste0(basePaths$rds, basePrefix, "_Boolean_Rules.rds"),
    boolnet      = paste0(basePaths$rds, basePrefix, "_boolnet.rds"),
    attractors   = paste0(basePaths$rds, basePrefix, "_attractors.rds"),
    attractorDf  = paste0(basePaths$rds, basePrefix, "_attractor_df.rds"),
    geneMap      = paste0(basePaths$rds, basePrefix, "_gene_map.rds"),
    
    # Attractor analysis files
    attractorEntropy = paste0(basePaths$rds, basePrefix, "_attractor_entropy.rds"),
    agingScore       = paste0(basePaths$rds, basePrefix, "_aging_score.rds"),
    
    # GraphML exports
    grnPlot                = paste0(basePaths$plots,   basePrefix, "_GRN.png"),
    grnPlotPreprocessed    = paste0(basePaths$plots,   basePrefix, "_GRN_preprocessed.png"),
    grnGraphml             = paste0(basePaths$graphml, basePrefix, "_GRN.graphml"),
    grnGraphmlPreprocessed = paste0(basePaths$graphml, basePrefix, "_GRN_preprocessed.graphml"),
    
    # Perturbation analysis (these would be generated in loops)
    perturbationResults = paste0(basePaths$rds, basePrefix, "_perturbation_results.rds"),
    rejuvenationTargets = paste0(basePaths$rds, basePrefix, "_rejuvenation_targets.rds")
  )
}

# =============================================================================
# Perturbation-Specific Paths
# =============================================================================

#' Build paths for individual gene perturbation results
#' @param basePaths List from buildBasePaths()
#' @param cellType Cell type string
#' @param trajectory Trajectory string
#' @param geneName Gene name
#' @param perturbationType "0" (knockout), "1" (overexpression), "KD_KD", "OE_OE", "KD_OE"
#' @return File path for perturbation result
getPerturbationFilePath <- function(basePaths, cellType, trajectory, geneName, perturbationType) {
  basePrefix <- paste0(cellType, "_", trajectory)
  paste0(basePaths$rds, basePrefix, "_", geneName, "_", perturbationType, ".rds")
}

#' Build paths for combination perturbation results
#' @param basePaths List from buildBasePaths()
#' @param cellType Cell type string
#' @param trajectory Trajectory string
#' @param gene1 First gene name
#' @param gene2 Second gene name
#' @param combinationType "KD_KD", "OE_OE", "KD_OE"
#' @return File path for combination perturbation result
getCombinationPerturbationFilePath <- function(basePaths, cellType, trajectory, gene1, gene2, combinationType) {
  basePrefix <- paste0(cellType, "_", trajectory)
  paste0(basePaths$rds, basePrefix, "_", gene1, "_", gene2, "_", combinationType, ".rds")
}

# =============================================================================
# Master Path Builder Function
# =============================================================================

#' Build comprehensive path structure for the project
#' @param config Configuration object from yaml::read_yaml("config.yaml")
#' @param cellType Optional cell type (if NULL, only static paths returned)
#' @param trajectory Optional trajectory (requires cellType)
#' @return Named list containing all relevant paths
buildProjectPaths <- function(config, cellType = NULL, trajectory = NULL) {
  basePaths <- buildBasePaths(config)
  
  paths <- list(
    base = basePaths,
    cellranger = buildCellrangerPaths(config),
    static = getStaticFilePaths(basePaths)
  )
  
  if (!is.null(cellType)) {
    paths$cellType <- getCellTypeFilePaths(basePaths, cellType)
    
    if (!is.null(trajectory)) {
      paths$trajectory <- getTrajectoryFilePaths(basePaths, cellType, trajectory)
    }
  }
  
  return(paths)
}

# =============================================================================
# Convenience Functions for Scripts
# =============================================================================

#' Initialize paths for a script with known cell type and trajectory
#' @param cellType Cell type string
#' @param trajectory Trajectory string
#' @return List of paths ready for use
initializeScriptPaths <- function(cellType = NULL, trajectory = NULL) {
  config <- yaml::read_yaml("config.yaml")
  buildProjectPaths(config, cellType, trajectory)
}

#' Initialize paths for scripts that need to select cell type/trajectory interactively
#' @param needsCellType Boolean indicating if cell type selection is needed
#' @param needsTrajectory Boolean indicating if trajectory selection is needed
#' @return List containing paths and selected cellType/trajectory
initializeInteractivePaths <- function(needsCellType = FALSE, needsTrajectory = FALSE) {
  config <- yaml::read_yaml("config.yaml")
  basePaths <- buildBasePaths(config)
  
  result <- list(
    paths      = buildProjectPaths(config),
    cellType   = NULL,
    trajectory = NULL
  )
  
  if (needsCellType) {
    if (!file.exists(paste0(basePaths$rds, "cell_types.rds"))) {
      stop("Cell types file not found. Run script 02 first.")
    }
    
    cellTypes       <- readRDS(paste0(basePaths$rds, "cell_types.rds"))
    result$cellType <- showCellTypeMenu(cellTypes)
    result$paths    <- buildProjectPaths(config, result$cellType)
    
    if (needsTrajectory) {
      trajectoryFile <- paste0(basePaths$rds, "retained_trajectories_", result$cellType, ".rds")
      if (!file.exists(trajectoryFile)) {
        stop("Trajectory file not found for ", result$cellType, ". Run script 03 first.")
      }
      
      trajectories      <- readRDS(trajectoryFile)
      result$trajectory <- showTrajectoryMenu(trajectories)
      result$paths      <- buildProjectPaths(config, result$cellType, result$trajectory)
    }
  }
  
  return(result)
}

# =============================================================================
# Directory Creation Helper
# =============================================================================

#' Ensure all required directories exist
#' @param paths Path structure from buildProjectPaths()
ensureProjectDirectories <- function(paths) {
  # Create all base directories
  lapply(paths$base, function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  })
  
  # Create monocle3 directories if they exist in the path structure
  if ("trajectory" %in% names(paths)) {
    monocle3Dirs <- paths$trajectory[grepl("monocle3", names(paths$trajectory))]
    lapply(monocle3Dirs, function(path) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    })
  }
  
  invisible(TRUE)
}
