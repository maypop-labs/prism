# =============================================================================
# pseudotimeManager.R
# Purpose: Functions for analyzing Monocle3 pseudotime trajectories
# =============================================================================

###############################################################################
# Core Analysis Functions
###############################################################################

#' Calculate Trajectory-Age Correlation
#'
#' Computes the Spearman correlation between pseudotime values and donor age
#' along a trajectory branch. Tests for monotonic relationship between cellular
#' progression (pseudotime) and chronological age.
#'
#' @param subCds A Monocle3 CellDataSet object containing a subset of cells
#'   from a single trajectory branch. Must have 'age' and 'donorID' columns
#'   in colData.
#'
#' @return A list containing:
#'   \item{correlation}{Spearman's rho correlation coefficient}
#'   \item{p_value}{P-value for the correlation test}
#'   \item{method}{Method used ("spearman_correlation" or "insufficient_data")}
#'   \item{n_cells}{Number of valid cells used in analysis}
#'   \item{n_donors}{Number of unique donors represented}
#'
#' @details
#' The function filters out cells with NA or non-finite pseudotime values.
#' Requires at least 10 valid cells to perform correlation analysis.
#' Uses Spearman correlation to test for monotonic (not necessarily linear)
#' relationship between pseudotime and age.
#'
#' @examples
#' \dontrun{
#' result <- calculateTrajectoryCorrelation(subCds)
#' if (result$correlation > 0.6 && result$p_value < 0.05) {
#'   message("Strong age-related trajectory detected")
#' }
#' }
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

#' Convert Seurat Object to Monocle3 CellDataSet
#'
#' Converts a Seurat v5 object to a Monocle3 CellDataSet format for trajectory
#' inference analysis. Extracts raw counts, cell metadata, and gene metadata.
#'
#' @param seuratObj A Seurat object (v5 format) containing single-cell
#'   RNA-seq data with counts layer and metadata.
#'
#' @return A Monocle3 CellDataSet object containing:
#'   \item{expression}{Raw count matrix}
#'   \item{colData}{Cell metadata from Seurat object}
#'   \item{rowData}{Gene metadata with gene_id and gene_short_name}
#'
#' @details
#' Extracts the "counts" layer from Seurat v5 LayerData structure.
#' Preserves all cell-level metadata from the original Seurat object.
#' Creates minimal gene metadata with gene identifiers.
#'
#' @examples
#' \dontrun{
#' seuratObj <- readRDS("celltype_seurat.rds")
#' cds <- convertSeuratToCDS(seuratObj)
#' }
convertSeuratToCDS <- function(seuratObj) {
  rawCounts <- LayerData(seuratObj, layer = "counts")
  genes <- rownames(rawCounts)
  cellMeta <- seuratObj@meta.data
  geneMeta <- data.frame(gene_id = genes, gene_short_name = genes, row.names = genes)
  new_cell_data_set(rawCounts, cellMeta, geneMeta)
}

#' Get Root Cells for Pseudotime Analysis
#'
#' Identifies cell barcodes corresponding to the youngest donor age for use
#' as root cells in trajectory inference. Root cells define the starting point
#' of pseudotime (pseudotime = 0).
#'
#' @param cds A Monocle3 CellDataSet object with 'age' column in colData.
#' @param rootAge Numeric value specifying the age to use as trajectory root.
#'   Typically the minimum age across all donors.
#'
#' @return Character vector of cell barcodes (rownames) for cells from donors
#'   matching the specified root age.
#'
#' @details
#' Root cells should represent the earliest/youngest state in the biological
#' process being modeled. Using youngest donor age ensures trajectories model
#' progression toward older cellular states.
#'
#' @examples
#' \dontrun{
#' rootAge <- min(config$ages)
#' rootCells <- getRootCells(cds, rootAge)
#' cds <- order_cells(cds, root_cells = rootCells)
#' }
getRootCells <- function(cds, rootAge) {
  meta <- as.data.frame(colData(cds))
  rownames(meta)[meta$age == rootAge]
}

#' Run Complete Pseudotime Analysis Pipeline
#'
#' Executes the full Monocle3 trajectory inference workflow: dimensionality
#' reduction, clustering, graph learning, and pseudotime ordering. Orders cells
#' along learned trajectories starting from specified root cells.
#'
#' @param cds A Monocle3 CellDataSet object with raw count data.
#' @param rootCells Character vector of cell barcodes to use as trajectory roots.
#'   Obtained from getRootCells().
#' @param verbose Logical indicating whether to print progress messages.
#'   Default is TRUE.
#'
#' @return A Monocle3 CellDataSet object with:
#'   \item{reducedDims}{PCA and UMAP embeddings}
#'   \item{clusters}{Leiden clustering assignments}
#'   \item{principal_graph}{Learned trajectory structure}
#'   \item{pseudotime}{Computed pseudotime values for each cell}
#'
#' @details
#' Pipeline steps:
#' \enumerate{
#'   \item PCA preprocessing (50 dimensions, log normalization with scaling)
#'   \item UMAP dimensionality reduction
#'   \item Leiden clustering
#'   \item Principal graph learning (trajectory structure)
#'   \item Cell ordering (pseudotime calculation from root cells)
#' }
#'
#' @examples
#' \dontrun{
#' rootCells <- getRootCells(cds, min(config$ages))
#' cds <- runPseudotime(cds, rootCells, verbose = TRUE)
#' pt <- pseudotime(cds)
#' }
runPseudotime <- function(cds, rootCells, verbose = TRUE, minBranchLen = 10) {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "log", scaling = TRUE, verbose = verbose)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose = verbose)
  cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden", verbose = verbose)
  cds <- learn_graph(cds, learn_graph_control = list(minimal_branch_len = minBranchLen))
  order_cells(cds, reduction_method = "UMAP", root_cells = rootCells)
}

###############################################################################
# Standardized Plotting Functions
###############################################################################

#' Create Standardized ggplot2 Plot
#'
#' Generates a ggplot2 object with standardized theming, geometry, and formatting
#' consistent across all PRISM visualizations. Provides centralized control over
#' plot aesthetics.
#'
#' @param plotData Data frame containing data to plot.
#' @param aes_mapping ggplot2 aesthetic mapping created with aes().
#' @param geom_type Character string specifying geometry type. Options:
#'   "point", "col", "point_colored". Default is "point".
#' @param title Character string for plot title. Default is empty.
#' @param subtitle Character string for plot subtitle. Default is empty.
#' @param config Configuration list containing plotting parameters:
#'   plotAlpha, pointSize, strokeSize, hjust.
#' @param additional_geoms List of additional ggplot2 geom layers to add.
#'   Default is NULL.
#' @param additional_themes Additional theme() elements to apply. Default is NULL.
#'
#' @return A ggplot2 object with standardized formatting.
#'
#' @details
#' Applies consistent styling:
#' \itemize{
#'   \item theme_minimal() base
#'   \item White background
#'   \item Centered titles (configurable via hjust)
#'   \item No minor grid lines
#'   \item Standardized point sizes and transparency
#' }
#'
#' @examples
#' \dontrun{
#' p <- createStandardPlot(
#'   plotData = df,
#'   aes_mapping = aes(x = pseudotime, y = age),
#'   geom_type = "point",
#'   title = "Age vs Pseudotime",
#'   config = config
#' )
#' }
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

#' Save Standardized Plot to File
#'
#' Saves a ggplot2 object using standardized dimensions and resolution settings
#' from configuration. Ensures consistent plot output across all PRISM analyses.
#'
#' @param plot A ggplot2 object to save.
#' @param filename Character string specifying output file path (including extension).
#' @param config Configuration list containing: figWidth, figHeight, figDPI.
#' @param width_scale Numeric multiplier for figure width. Default is 1.
#' @param height_scale Numeric multiplier for figure height. Default is 1.
#'
#' @return NULL (invisibly). Writes plot file to disk.
#'
#' @details
#' Uses ggsave() with configuration-defined defaults:
#' \itemize{
#'   \item Width and height in inches (scalable via multipliers)
#'   \item DPI resolution from config
#'   \item Output format determined by filename extension
#' }
#'
#' @examples
#' \dontrun{
#' saveStandardPlot(
#'   plot = p,
#'   filename = "trajectory_plot.png",
#'   config = config,
#'   width_scale = 1.2,
#'   height_scale = 1.0
#' )
#' }
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

#' Plot Age versus Pseudotime for Single Trajectory
#'
#' Creates a scatter plot showing the relationship between donor age and
#' pseudotime for cells in a single trajectory branch. Includes linear fit,
#' correlation statistics, and donor-level coloring.
#'
#' @param subCds A Monocle3 CellDataSet containing cells from one trajectory branch.
#' @param leafName Character string identifying the trajectory leaf node.
#' @param correlationResult List containing correlation statistics from
#'   calculateTrajectoryCorrelation(): correlation, p_value, n_cells, n_donors.
#' @param config Configuration list with plotting parameters.
#'
#' @return A ggplot2 object showing age vs pseudotime scatter plot with:
#'   \item{points}{Colored by donor ID}
#'   \item{fit line}{Linear regression with confidence interval}
#'   \item{subtitle}{Correlation coefficient, p-value, cell count, donor count}
#'
#' @details
#' Uses viridis "turbo" color palette for donor visualization.
#' Removes legend as donor colors are for visual grouping only.
#' Filters out cells with invalid pseudotime values.
#'
#' @examples
#' \dontrun{
#' corResult <- calculateTrajectoryCorrelation(subCds)
#' p <- plotAgeVsPseudotime(subCds, "Y_123", corResult, config)
#' }
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

#' Create Plots for All Retained Trajectories
#'
#' Generates age vs pseudotime plots for all trajectory branches that passed
#' validation criteria. Creates both individual plots per trajectory and a
#' combined multi-panel figure when multiple trajectories exist.
#'
#' @param cds Full Monocle3 CellDataSet object with trajectory graph.
#' @param retainedLeaves Character vector of leaf node names for valid trajectories.
#' @param branchStats Data frame with trajectory statistics (leaf, corWithAge,
#'   pValue, nCells, nDonors).
#' @param basePaths List of base directory paths from pathManager.
#' @param cellType Character string specifying cell type being analyzed.
#' @param config Configuration list with plotting and analysis parameters.
#' @param saveIndividual Logical indicating whether to save individual trajectory
#'   plots. Default is TRUE.
#'
#' @return Named list of ggplot2 objects, one per retained trajectory.
#'   Names are leaf node identifiers.
#'
#' @details
#' For each retained trajectory:
#' \enumerate{
#'   \item Extracts cells belonging to trajectory path from root to leaf
#'   \item Creates age vs pseudotime scatter plot with statistics
#'   \item Saves individual plot (if saveIndividual = TRUE)
#' }
#'
#' If multiple trajectories exist:
#' \itemize{
#'   \item Arranges plots in grid (max 3 columns)
#'   \item Adds overall title with cell type and trajectory count
#'   \item Saves combined multi-panel figure
#' }
#'
#' @examples
#' \dontrun{
#' plots <- createTrajectoryPlots(
#'   cds = cds,
#'   retainedLeaves = c("Y_123", "Y_456"),
#'   branchStats = branch_stats,
#'   basePaths = paths$base,
#'   cellType = "Fibroblasts",
#'   config = config,
#'   saveIndividual = TRUE
#' )
#' }
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
