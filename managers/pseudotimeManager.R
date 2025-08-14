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