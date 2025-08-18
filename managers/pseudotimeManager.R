# =============================================================================
# pseudotimeManager.R
# Purpose: Functions for analyzing Monocle3 pseudotime trajectories
# =============================================================================

###############################################################################
# Core Analysis Functions
###############################################################################

calculateTrajectoryCorrelation <- function(subCds) {
  # Extract data
  pt <- pseudotime(subCds)
  ages <- colData(subCds)$age
  donors <- colData(subCds)$donorID
  
  # Clean data
  valid_idx <- !is.na(pt) & is.finite(pt) & !is.na(ages)
  
  if (sum(valid_idx) < 10) {
    return(list(
      correlation = NA, 
      p_value = NA, 
      method = "insufficient_data",
      n_cells = 0,
      n_donors = 0
    ))
  }
  
  pt_clean <- pt[valid_idx]
  ages_clean <- ages[valid_idx]
  donors_clean <- donors[valid_idx]
  n_donors <- length(unique(donors_clean))
  
  # Spearman correlation test
  spearman_test <- cor.test(pt_clean, ages_clean, method = "spearman")
  
  return(list(
    correlation = spearman_test$estimate,
    p_value = spearman_test$p.value,
    method = "spearman_correlation",
    n_cells = sum(valid_idx),
    n_donors = n_donors
  ))
}

convertSeuratToCDS <- function(seuratObj) {
  rawCounts <- LayerData(seuratObj, layer = "counts")
  genes <- rownames(rawCounts)
  cellMeta <- seuratObj@meta.data
  geneMeta <- data.frame(gene_id = genes, gene_short_name = genes, row.names = genes)
  new_cell_data_set(rawCounts, cellMeta, geneMeta)
}

getRootCells <- function(cds, rootAge) {
  meta <- as.data.frame(colData(cds))
  rownames(meta)[meta$age == rootAge]
}

runPseudotime <- function(cds, rootCells, verbose = TRUE) {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "log", scaling = TRUE, verbose = verbose)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose = verbose)
  cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden", verbose = verbose)
  cds <- learn_graph(cds)
  order_cells(cds, reduction_method = "UMAP", root_cells = rootCells)
}

# =============================================================================
# smoothPseudotimeExpression
# 
# Smooths gene expression profiles along pseudotime using a moving average window
# and saves the updated Monocle3 object with smoothed expression as a new assay.
# =============================================================================

#' Smooth Gene Expression Along Pseudotime
#'
#' Applies moving average smoothing to gene expression profiles along pseudotime
#' to reduce single-cell noise while preserving temporal trends. The smoothed
#' expression is stored as a new assay called "smoothed_expr" in the Monocle3 object.
#'
#' @param cellType Character string specifying the cell type to process.
#' @param trajectory Character string specifying the trajectory name (e.g., "Branch_1").
#' @param winSizePercent Numeric value between 0 and 1 specifying the smoothing 
#'   window size as a percentage of total cells. Default is 0.01 (1%).
#' @param saveResults Logical indicating whether to save the smoothed Monocle3 object.
#'   Default is TRUE.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param basePath Character string specifying the base project path. If NULL (default),
#'   uses interactive path initialization.
#'
#' @return A Monocle3 CellDataSet object with smoothed expression added as the
#'   "smoothed_expr" assay.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads the Monocle3 object from the specified trajectory directory
#'   \item Computes pseudotime if not already present
#'   \item Orders cells by pseudotime and applies moving average smoothing
#'   \item Stores smoothed values as a new assay called "smoothed_expr"
#'   \item Saves the updated object to the monocle3Smoothed directory
#' }
#'
#' The moving average window size is calculated as max(5, floor(winSizePercent * nCells))
#' to ensure a minimum window of 5 cells while scaling with dataset size.
#'
#' @examples
#' \dontrun{
#' # Basic usage with interactive path selection
#' cds <- smoothPseudotimeExpression(cellType = "Keratinocyte", 
#'                                   trajectory = "Branch_1")
#'
#' # Custom window size (15% of cells)
#' cds <- smoothPseudotimeExpression(cellType = "Keratinocyte",
#'                                   trajectory = "Branch_1", 
#'                                   winSizePercent = 0.15)
#'
#' # Specify custom base path
#' cds <- smoothPseudotimeExpression(cellType = "Keratinocyte",
#'                                   trajectory = "Branch_1",
#'                                   basePath = "/path/to/project")
#' }
#'
#' @seealso
#' \code{\link{load_monocle_objects}}, \code{\link{save_monocle_objects}},
#' \code{\link{getCellTypeFilePaths}}, \code{\link{getTrajectoryFilePaths}}
#'
#' @export
smoothPseudotimeExpression <- function(cellType, 
                                       trajectory, 
                                       winSizePercent = 0.10, 
                                       saveResults = TRUE, 
                                       verbose = TRUE,
                                       basePath = NULL) {
  
  # --- Parameter Validation ---
  if (missing(cellType) || missing(trajectory)) {
    stop("Both cellType and trajectory must be specified")
  }
  
  if (winSizePercent <= 0 || winSizePercent > 1) {
    stop("winSizePercent must be between 0 and 1")
  }
  
  # --- Initialize Paths ---
  if (is.null(basePath)) {
    pathInfo <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
    paths <- pathInfo$paths
  } else {
    paths <- list(base = basePath)
  }

  ptPaths <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
  
  # --- Validate Input Directory ---
  if (!dir.exists(ptPaths$monocle3)) {
    stop("Monocle3 directory not found: ", ptPaths$monocle3)
  }
  
  # --- Load Monocle3 Object ---
  if (verbose) {
    message("Loading Monocle3 object from: ", ptPaths$monocle3)
  }
  
  cds <- load_monocle_objects(directory_path = ptPaths$monocle3)
  
  # --- Ensure Pseudotime is Computed ---
  if (is.null(colData(cds)$Pseudotime)) {
    if (verbose) {
      message("Computing pseudotime")
    }
    colData(cds)$Pseudotime <- pseudotime(cds)
  }
  
  # --- Moving Average Smoothing ---
  cellOrder <- order(colData(cds)$Pseudotime)
  nCells    <- length(cellOrder)
  winSize   <- max(5, floor(winSizePercent * nCells))
  halfWin   <- floor(winSize / 2)
  exprMat   <- assay(cds, "counts")
  
  # Initialize smoothed matrix with same dimensions and names
  smoothMat <- matrix(0, 
                      nrow = nrow(exprMat), 
                      ncol = ncol(exprMat), 
                      dimnames = dimnames(exprMat))
  
  if (verbose) {
    message("Smoothing expression values with window size: ", winSize, 
            " (", round(winSizePercent * 100, 1), "% of ", nCells, " cells)")
  }
  
  # Apply moving average smoothing
  for (i in seq_len(nCells)) {
    startIdx <- max(1, i - halfWin)
    endIdx   <- min(nCells, i + halfWin)
    idx      <- cellOrder[startIdx:endIdx]
    smoothMat[, cellOrder[i]] <- rowMeans(exprMat[, idx, drop = FALSE])
  }
  
  # Add smoothed expression as new assay
  assay(cds, "smoothed_expr") <- smoothMat
  
  # --- Save Results ---
  if (saveResults) {
    if (verbose) {
      message("Saving smoothed Monocle3 object to: ", ptPaths$monocle3Smoothed)
    }
    save_monocle_objects(cds = cds, 
                         directory_path = ptPaths$monocle3Smoothed, 
                         comment = cellType)
  }
  
  if (verbose) {
    message("Smoothing completed successfully!")
  }
  
  # Return the updated CDS object
  return(cds)
}

###############################################################################
# Standardized Plotting Functions (Option A)
###############################################################################

createStandardPlot <- function(plotData, aes_mapping, geom_type = "point", 
                               title = "", subtitle = "", config, 
                               additional_geoms = NULL, additional_themes = NULL) {
  p <- ggplot(plotData, aes_mapping)
  
  # Add geometry based on type
  if (geom_type == "point") {
    p <- p + geom_point(
      alpha = config$plotAlpha,
      size = config$pointSize,
      stroke = config$strokeSize
    )
  } else if (geom_type == "col") {
    p <- p + geom_col(alpha = config$plotAlpha)
  } else if (geom_type == "point_colored") {
    p <- p + geom_point(
      alpha = config$plotAlpha,
      size = config$pointSize,
      stroke = config$strokeSize
    )
  }
  
  # Add additional geoms if provided
  if (!is.null(additional_geoms)) {
    for (geom in additional_geoms) {
      p <- p + geom
    }
  }
  
  # Standard theme and formatting
  p <- p + 
    labs(title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = config$hjust, size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Add additional theme elements if provided
  if (!is.null(additional_themes)) {
    p <- p + additional_themes
  }
  
  return(p)
}

saveStandardPlot <- function(plot, filename, config, width_scale = 1, height_scale = 1) {
  ggsave(
    filename = filename,
    plot = plot,
    width = config$figWidth * width_scale,
    height = config$figHeight * height_scale,
    dpi = config$figDPI,
    units = "in"
  )
}

###############################################################################
# Visualization Functions
###############################################################################

# Function to create scatter plot for a single trajectory
plotAgeVsPseudotime <- function(subCds, leafName, correlationResult, config) {
  # Extract data
  pt       <- pseudotime(subCds)
  ages     <- colData(subCds)$age
  donorIDs <- colData(subCds)$donorID
  
  # Create data frame
  plotData <- data.frame(
    pseudotime = pt,
    age        = ages,
    donorID    = donorIDs
  ) %>%
    filter(!is.na(pseudotime) & is.finite(pseudotime) & !is.na(age))
  
  # Define additional geoms for this specific plot
  additional_geoms <- list(
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1),
    scale_color_viridis_d(option = "turbo")
  )
  
  # Define additional theme elements
  additional_themes <- theme(
    legend.position = "none"
  )
  
  # Create the plot using standardized function
  p <- createStandardPlot(
    plotData = plotData,
    aes_mapping = aes(x = pseudotime, y = age, color = as.factor(donorID)),
    geom_type = "point_colored",
    title = paste0("Trajectory ", leafName),
    subtitle = paste0("r = ", round(correlationResult$correlation, 3), 
                      ", p = ", format(correlationResult$p_value, scientific = TRUE, digits = 2),
                      "\nCells: ", correlationResult$n_cells, 
                      ", Donors: ", correlationResult$n_donors),
    config = config,
    additional_geoms = additional_geoms,
    additional_themes = additional_themes
  ) +
    labs(
      x = "Pseudotime",
      y = "Donor Age",
      color = "donorID"
    )
  
  return(p)
}

# Function to create plots for all retained trajectories
createTrajectoryPlots <- function(cds, retainedLeaves, branchStats, 
                                  basePaths, cellType, config, saveIndividual = TRUE) {
  
  # Get trajectory graph info
  trajGraph <- principal_graph(cds)[["UMAP"]]
  closestVertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  rootCells <- getRootCells(cds, min(colData(cds)$age))
  rootVertices <- unique(V(trajGraph)$name[closestVertex[rootCells, 1]])
  
  plotList <- list()
  
  # Create plot for each retained trajectory
  for (leaf in retainedLeaves) {
    message("Creating plot for trajectory ", leaf)
    ptPaths <- getTrajectoryFilePaths(basePaths, cellType, leaf)
    
    # Extract cells for this trajectory
    sp <- shortest_paths(trajGraph, from = rootVertices, to = leaf, weights = NA)$vpath[[1]]
    pathNodeIdx <- as.integer(sub("Y_", "", as_ids(sp)))
    cellBarcodes <- rownames(closestVertex)[closestVertex[, 1] %in% pathNodeIdx]
    subCds <- cds[, cellBarcodes]
    
    # Get correlation results for this trajectory
    leafStats <- branchStats[branchStats$leaf == leaf, ]
    correlationResult <- list(
      correlation = leafStats$corWithAge,
      p_value     = leafStats$pValue,
      n_cells     = leafStats$nCells,
      n_donors    = leafStats$nDonors
    )
    
    # Create plot using standardized function
    p <- plotAgeVsPseudotime(subCds, leaf, correlationResult, config)
    plotList[[leaf]] <- p
    
    # Save individual plot if requested
    if (saveIndividual) {
      saveStandardPlot(
        plot = p,
        filename = ptPaths$ageVsPseudotimePlot,
        config = config,
        width_scale = 1.2,
        height_scale = 1.0
      )
    }
  }
  
  # Create combined plot if multiple trajectories
  if (length(plotList) > 1) {
    # Arrange plots in a grid
    nPlots <- length(plotList)
    nCols <- min(3, nPlots)  # Max 3 columns
    nRows <- ceiling(nPlots / nCols)
    
    combinedPlot <- wrap_plots(plotList, ncol = nCols)
    
    # Add overall title
    combinedPlot <- combinedPlot + 
      plot_annotation(
        title = paste0("Age vs Pseudotime: ", cellType, " (", nPlots, " Retained Trajectories)"),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = config$hjust))
      )
    
    # Save combined plot using standardized function
    saveStandardPlot(
      plot = combinedPlot,
      filename = ptPaths$ageVsPseudotimeCombinedPlot,
      config = config,
      width_scale = nCols * 0.92,
      height_scale = nRows * 0.77
    )
    
    message("Saved combined plot with ", nPlots, " trajectories")
  }
  
  return(plotList)
}