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
#' @export
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
    graphml    = paste0(rootPath, "results/graphml/"),
    scenic     = paste0(rootPath, "results/scenic/")
  )
}

#' Build CellRanger input paths
#' @param config Configuration object
#' @return Vector of CellRanger paths for each sample
#' @export
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
#' @export
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
#' @export
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
#' @export
getTrajectoryFilePaths <- function(basePaths, cellType, trajectory) {
  basePrefix <- paste0(cellType, "_", trajectory)
  
  list(
    # Monocle3 directories
    monocle3             = paste0(basePaths$monocle3, "monocle3_", basePrefix, "/"),
    monocle3GeneSwitches = paste0(basePaths$monocle3, "monocle3_", basePrefix, "_geneSwitches/"),

    # SCENIC directory
    scenic = paste0(basePaths$scenic, basePrefix, "/"),

    # Pseudotime Analysis
    ageVsPseudotimePlot         = paste0(basePaths$plots, basePrefix, "_age_vs_pseudotime.png"),
    ageVsPseudotimeCombinedPlot = paste0(basePaths$plots, basePrefix, "_age_vs_pseudotime_combined.png"),
    
    # Analysis intermediate files
    geneMap         = paste0(basePaths$rds, basePrefix, "_gene_mapping.rds"),
    geneSwitches    = paste0(basePaths$rds, basePrefix, "_geneSwitches.rds"),
    geneSwitchesTsv = paste0(basePaths$tsv, basePrefix, "_geneSwitches.tsv"),
    grn             = paste0(basePaths$rds, basePrefix, "_GRN.rds"),
    grnEdges        = paste0(basePaths$rds, basePrefix, "_GRN_edges.rds"),
    grnEdgesTsv     = paste0(basePaths$tsv, basePrefix, "_GRN_edges.tsv"),

    # Boolean network files
    booleanRules = paste0(basePaths$rds, basePrefix, "_Boolean_Rules.rds"),
    boolNet      = paste0(basePaths$rds, basePrefix, "_boolNet.rds"),
    attractors   = paste0(basePaths$rds, basePrefix, "_attractors.rds"),
    attractorDf  = paste0(basePaths$rds, basePrefix, "_attractor_df.rds"),
    geneMap      = paste0(basePaths$rds, basePrefix, "_gene_map.rds"),
    boolNetTsv   = paste0(basePaths$tsv, basePrefix, "_boolNet_analysis.tsv"),
    
    # Attractor analysis files
    attractorEntropy        = paste0(basePaths$rds, basePrefix, "_attractor_entropy.rds"),
    agingScore              = paste0(basePaths$rds, basePrefix, "_aging_score.rds"),
    referenceVectors        = paste0(basePaths$rds, basePrefix, "_reference_vectors.rds"),
    attractorScoresCombined = paste0(basePaths$rds, basePrefix, "_attractor_scores_combined.rds"),
    networkAgingSummary     = paste0(basePaths$rds, basePrefix, "_network_aging_summary.rds"),

    # Perturbation analysis files
    singleTargetsKD         = paste0(basePaths$rds, basePrefix, "_single_targets_KD.rds"),
    singleTargetsOE         = paste0(basePaths$rds, basePrefix, "_single_targets_OE.rds"),
    perturbationSummary     = paste0(basePaths$rds, basePrefix, "_perturbation_summary.rds"),

    # Double perturbation analysis files
    doubleTargetsKD         = paste0(basePaths$rds, basePrefix, "_double_targets_KD_KD.rds"),
    doubleTargetsOE         = paste0(basePaths$rds, basePrefix, "_double_targets_OE_OE.rds"),
    doubleTargetsMix        = paste0(basePaths$rds, basePrefix, "_double_targets_KD_OE.rds"),

    # Final analysis files
    finalTargetRankings     = paste0(basePaths$rds, basePrefix, "_final_target_rankings.rds"),
    targetSummaryTsv        = paste0(basePaths$tsv, basePrefix, "_target_summary.tsv"),

    # Boolean analysis reports (additional TSV outputs)
    booleanRulesAnalysis    = paste0(basePaths$tsv, basePrefix, "_boolean_rules_analysis.tsv"),
    booleanRulesFlat        = paste0(basePaths$tsv, basePrefix, "_boolean_rules.tsv"),
    
    # GraphML exports
    grnPlot    = paste0(basePaths$plots,   basePrefix, "_GRN.png"),
    grnGraphml = paste0(basePaths$graphml, basePrefix, "_GRN.graphml"),

    # Perturbation analysis (these would be generated in loops)
    perturbationResults = paste0(basePaths$rds, basePrefix, "_perturbation_results.rds"),
    rejuvenationTargets = paste0(basePaths$rds, basePrefix, "_rejuvenation_targets.rds")
  )
}

# =============================================================================
# Master Path Builder Function
# =============================================================================

#' Build comprehensive path structure for the project
#' @param config Configuration object from yaml::read_yaml("config.yaml")
#' @param cellType Optional cell type (if NULL, only static paths returned)
#' @param trajectory Optional trajectory (requires cellType)
#' @return Named list containing all relevant paths
#' @export
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

#' Initialize paths for scripts that need to select cell type/trajectory interactively
#' @param needsCellType Boolean indicating if cell type selection is needed
#' @param needsTrajectory Boolean indicating if trajectory selection is needed
#' @return List containing paths and selected cellType/trajectory
#' @export
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
#' @export
ensureProjectDirectories <- function(paths) {
  # Create all base directories
  lapply(paths$base, function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  })
  
  # Create trajectory-specific directories if they exist in the path structure
  if ("trajectory" %in% names(paths)) {
    trajectoryDirs <- paths$trajectory[grepl("monocle3|scenic", names(paths$trajectory))]
    lapply(trajectoryDirs, function(path) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    })
  }
  
  invisible(TRUE)
}

# =============================================================================
# Universal File I/O Functions
# =============================================================================

#' Universal save function with automatic file type detection
#' @param object Object to save
#' @param filepath Full file path with extension
#' @param config Configuration object
#' @param objectName Optional name for verbose output (defaults to deparse)
#' @export
saveObject <- function(object, filepath, config, objectName = NULL) {
  if (is.null(objectName)) {
    objectName <- deparse(substitute(object))
  }
  
  # Auto-detect file type from extension
  ext <- tools::file_ext(filepath)
  
  if (config$verbose) {
    message("Saving ", objectName, " (", toupper(ext), ")")
  }
  
  switch(tolower(ext),
    "rds" = saveRDS(object, file = filepath),
    "tsv" = write.table(object, file = filepath, sep = "\t", quote = FALSE, row.names = FALSE),
    "txt" = writeLines(object, con = filepath),
    "json" = jsonlite::write_json(object, path = filepath),
    "graphml" = igraph::write_graph(object, file = filepath, format = "graphml"),
    stop("Unsupported file type: ", ext)
  )
}

#' Universal load function with automatic file type detection  
#' @param filepath Full file path with extension
#' @param config Configuration object
#' @param objectName Optional name for verbose output
#' @param required Whether to stop on missing file (default TRUE)
#' @return Loaded object
#' @export
loadObject <- function(filepath, config, objectName = NULL, required = TRUE) {
  if (!file.exists(filepath)) {
    msg <- paste("File not found:", filepath)
    if (required) stop(msg) else {warning(msg); return(NULL)}
  }
  
  if (is.null(objectName)) {
    objectName <- tools::file_path_sans_ext(basename(filepath))
  }
  
  ext <- tools::file_ext(filepath)
  
  if (config$verbose) {
    message("Loading ", objectName, " (", toupper(ext), ")")
  }
  
  switch(tolower(ext),
    "rds" = readRDS(filepath),
    "tsv" = read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    "txt" = readLines(filepath),
    "json" = jsonlite::read_json(filepath),
    "graphml" = igraph::read_graph(filepath, format = "graphml"),
    stop("Unsupported file type: ", ext)
  )
}

# =============================================================================
# Special Case Handlers for Monocle3
# =============================================================================

#' Load Monocle3 objects (special case for directories)
#' @param directory Directory path containing Monocle3 objects
#' @param config Configuration object
#' @param objectName Optional name for verbose output
#' @return Monocle3 CellDataSet object
#' @export
loadMonocle3 <- function(directory, config, objectName = "Monocle3 object") {
  if (!dir.exists(directory)) stop("Directory not found: ", directory)
  if (config$verbose) message("Loading ", objectName)
  load_monocle_objects(directory_path = directory)
}

#' Save Monocle3 objects (special case for directories)
#' @param cds Monocle3 CellDataSet object
#' @param directory Directory path to save Monocle3 objects
#' @param config Configuration object
#' @param objectName Optional name for verbose output
#' @export
saveMonocle3 <- function(cds, directory, config, objectName = "Monocle3 object") {
  if (config$verbose) message("Saving ", objectName)
  save_monocle_objects(cds = cds, directory_path = directory)
}

# =============================================================================
# Convenience Wrapper Functions (for backward compatibility)
# =============================================================================

# Boolean Rules
#' Load Boolean rules from file
#' @param ptPaths Path structure containing booleanRules path
#' @param config Configuration object
#' @return Boolean rules object
#' @export
loadBooleanRules <- function(ptPaths, config) {
  loadObject(ptPaths$booleanRules, config, "Boolean rules")
}

#' Save Boolean rules to file
#' @param boolRules Boolean rules object to save
#' @param ptPaths Path structure containing booleanRules path
#' @param config Configuration object
#' @export
saveBooleanRules <- function(boolRules, ptPaths, config) {
  saveObject(boolRules, ptPaths$booleanRules, config, "Boolean rules")
}

# GRN Edges
#' Load GRN edges from file
#' @param ptPaths Path structure containing grnEdges path
#' @param config Configuration object
#' @return GRN edges data frame
#' @export
loadGrnEdges <- function(ptPaths, config) {
  loadObject(ptPaths$grnEdges, config, "GRN edges")
}

# Gene Mapping
#' Save gene mapping to file
#' @param geneMap Gene mapping object to save
#' @param ptPaths Path structure containing geneMap path
#' @param config Configuration object
#' @export
saveGeneMap <- function(geneMap, ptPaths, config) {
  saveObject(geneMap, ptPaths$geneMap, config, "gene mapping")
}

# Gene Switches
#' Load switch genes from file
#' @param ptPaths Path structure containing geneSwitches path
#' @param config Configuration object
#' @return Switch genes data frame
#' @export
loadSwitchGenes <- function(ptPaths, config) {
  loadObject(ptPaths$geneSwitches, config, "switch genes")
}

#' Save switch genes to file
#' @param switchGenes Switch genes object to save
#' @param ptPaths Path structure containing geneSwitches path
#' @param config Configuration object
#' @export
saveGeneSwitches <- function(switchGenes, ptPaths, config) {
  saveObject(switchGenes, ptPaths$geneSwitches, config, "switch genes")
}

# Seurat Objects
#' Load merged Seurat object from file
#' @param paths Path structure containing seuratMerged path
#' @param config Configuration object
#' @return Merged Seurat object
#' @export
loadMergedSeurat <- function(paths, config) {
  loadObject(paths$static$seuratMerged, config, "merged Seurat object")
}

#' Save merged Seurat object to file
#' @param seuratMerged Merged Seurat object to save
#' @param paths Path structure containing seuratMerged path
#' @param config Configuration object
#' @export
saveMergedSeurat <- function(seuratMerged, paths, config) {
  saveObject(seuratMerged, paths$static$seuratMerged, config, "merged Seurat object")
}

#' Load cell type-specific Seurat object from file
#' @param ctPaths Cell type path structure containing seuratObject path
#' @param config Configuration object
#' @return Cell type Seurat object
#' @export
loadSeuratByCellType <- function(ctPaths, config) {
  loadObject(ctPaths$seuratObject, config, "cell type Seurat object")
}

#' Save cell type-specific Seurat object to file
#' @param seuratObject Seurat object to save
#' @param ctPaths Cell type path structure containing seuratObject path
#' @param config Configuration object
#' @export
saveSeuratByCellType <- function(seuratObject, ctPaths, config) {
  saveObject(seuratObject, ctPaths$seuratObject, config, "cell type Seurat object")
}

# Monocle3 Objects
#' Load pseudotime trajectory from Monocle3 directory
#' @param ptPaths Path structure containing monocle3 path
#' @param config Configuration object
#' @return Monocle3 CellDataSet object
#' @export
loadPseudotimeTrajectory <- function(ptPaths, config) {
  loadMonocle3(ptPaths$monocle3, config, "pseudotime trajectory")
}

#' Save pseudotime trajectory to Monocle3 directory
#' @param cds Monocle3 CellDataSet object to save
#' @param ptPaths Path structure containing monocle3 path
#' @param config Configuration object
#' @export
savePseudotimeTrajectory <- function(cds, ptPaths, config) {
  saveMonocle3(cds, ptPaths$monocle3, config, "pseudotime trajectory")
}

#' Load GeneSwitches trajectory from Monocle3 directory
#' @param ptPaths Path structure containing monocle3GeneSwitches path
#' @param config Configuration object
#' @return Monocle3 CellDataSet object
#' @export
loadMonocle3GeneSwitches <- function(ptPaths, config) {
  loadMonocle3(ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
}

#' Save GeneSwitches trajectory to Monocle3 directory
#' @param cds Monocle3 CellDataSet object to save
#' @param ptPaths Path structure containing monocle3GeneSwitches path
#' @param config Configuration object
#' @export
saveMonocle3GeneSwitches <- function(cds, ptPaths, config) {
  saveMonocle3(cds, ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
}

# Reports and Analysis
#' Save GeneSwitches report to file
#' @param switchOut GeneSwitches output object to save
#' @param ptPaths Path structure containing geneSwitchesTsv path
#' @param config Configuration object
#' @export
saveGeneSwitchesReport <- function(switchOut, ptPaths, config) {
  saveObject(switchOut, ptPaths$geneSwitchesTsv, config, "gene switches report")
}

#' Save BoolNet analysis report to file
#' @param boolNetReport BoolNet report object to save
#' @param ptPaths Path structure containing boolNetTsv path
#' @param config Configuration object
#' @export
saveBoolNetReport <- function(boolNetReport, ptPaths, config) {
  saveObject(boolNetReport, ptPaths$boolNetTsv, config, "BoolNet analysis report")
}

# BoolNet Objects
#' Save BoolNet network to file
#' @param boolNetwork BoolNet network object to save
#' @param ptPaths Path structure containing boolNet path
#' @param config Configuration object
#' @export
saveBoolNetwork <- function(boolNetwork, ptPaths, config) {
  saveObject(boolNetwork, ptPaths$boolNet, config, "BoolNet network")
}

#' Save attractors to file
#' @param attractors Attractors object to save
#' @param ptPaths Path structure containing attractors path
#' @param config Configuration object
#' @export
saveAttractors <- function(attractors, ptPaths, config) {
  saveObject(attractors, ptPaths$attractors, config, "attractors")
}
