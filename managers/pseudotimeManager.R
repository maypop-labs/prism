# =============================================================================
# pseudotimeManager.R
# Purpose: Functions for analyzing Monocle3 pseudotime trajectories
# =============================================================================

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

runPseudotime <- function(cds, rootCells) {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "log", scaling = TRUE, verbose = config$verbose)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose = config$verbose)
  cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden", verbose = config$verbose)
  cds <- learn_graph(cds)
  order_cells(cds, reduction_method = "UMAP", root_cells = rootCells)
}

###############################################################################
# Visualization Functions
###############################################################################

# Function to create scatter plot for a single trajectory
plotAgeVsPseudotime <- function(subCds, leafName, correlationResult) {
  # Extract data
  pt <- pseudotime(subCds)
  ages <- colData(subCds)$age
  donorIDs <- colData(subCds)$donorID
  
  # Create data frame
  plotData <- data.frame(
    pseudotime = pt,
    age = ages,
    donorID = donorIDs
  ) %>%
    filter(!is.na(pseudotime) & is.finite(pseudotime) & !is.na(age))
  
  # Create the plot
  p <- ggplot(plotData, aes(x = pseudotime, y = age)) +
    geom_point(aes(color = as.factor(donorID)), alpha = 0.6, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    labs(
      title = paste0("Trajectory ", leafName),
      subtitle = paste0("r = ", round(correlationResult$correlation, 3), 
                        ", p = ", format(correlationResult$p_value, scientific = TRUE, digits = 2),
                        "\nCells: ", correlationResult$n_cells, 
                        ", Donors: ", correlationResult$n_donors),
      x = "Pseudotime",
      y = "Donor Age",
      color = "donorID"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "none",  # Remove legend to save space
      panel.grid.minor = element_blank()
    ) +
    scale_color_viridis_d(option = "turbo")
  
  return(p)
}

# Function to create plots for all retained trajectories
createTrajectoryPlots <- function(cds, retainedLeaves, branchStats, 
                                  plotPath, cellType, saveIndividual = TRUE) {
  
  # Get trajectory graph info
  trajGraph <- principal_graph(cds)[["UMAP"]]
  closestVertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  rootCells <- getRootCells(cds, min(colData(cds)$age))
  rootVertices <- unique(V(trajGraph)$name[closestVertex[rootCells, 1]])
  
  plotList <- list()
  
  # Create plot for each retained trajectory
  for (leaf in retainedLeaves) {
    message("Creating plot for trajectory ", leaf)
    
    # Extract cells for this trajectory
    sp <- shortest_paths(trajGraph, from = rootVertices, to = leaf, weights = NA)$vpath[[1]]
    pathNodeIdx <- as.integer(sub("Y_", "", as_ids(sp)))
    cellBarcodes <- rownames(closestVertex)[closestVertex[, 1] %in% pathNodeIdx]
    subCds <- cds[, cellBarcodes]
    
    # Get correlation results for this trajectory
    leafStats <- branchStats[branchStats$leaf == leaf, ]
    correlationResult <- list(
      correlation = leafStats$corWithAge,
      p_value = leafStats$pValue,
      n_cells = leafStats$nCells,
      n_donors = leafStats$nDonors
    )
    
    # Create plot
    p <- plotAgeVsPseudotime(subCds, leaf, correlationResult)
    plotList[[leaf]] <- p
    
    # Save individual plot if requested
    if (saveIndividual) {
      ggsave(
        filename = paste0(plotPath, "age_vs_pseudotime_", cellType, "_", leaf, ".png"),
        plot = p,
        width = 8, height = 6, dpi = 300
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
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    # Save combined plot
    ggsave(
      filename = paste0(plotPath, "age_vs_pseudotime_", cellType, "_combined.png"),
      plot = combinedPlot,
      width = nCols * 6, height = nRows * 5, dpi = 300
    )
    
    message("Saved combined plot with ", nPlots, " trajectories")
  }
  
  return(plotList)
}

# Function to create summary statistics plot
createSummaryPlot <- function(branchStats, plotPath, cellType) {
  # Filter to retained trajectories (assuming correlation > threshold)
  retainedStats <- branchStats %>%
    filter(!is.na(corWithAge) & corWithAge >= 0.6) %>%  # Adjust threshold as needed
    arrange(desc(corWithAge))
  
  if (nrow(retainedStats) == 0) {
    message("No retained trajectories to plot")
    return(NULL)
  }
  
  # Create summary bar plot
  summaryPlot <- ggplot(retainedStats, aes(x = reorder(leaf, corWithAge), y = corWithAge)) +
    geom_col(aes(fill = corWithAge), alpha = 0.8) +
    geom_text(aes(label = paste0("n=", nCells)), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = paste0("Trajectory Age Correlations: ", cellType),
      subtitle = paste0(nrow(retainedStats), " retained trajectories"),
      x = "Trajectory",
      y = "Age-Pseudotime Correlation",
      fill = "Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    ) +
    scale_fill_gradient2(low = "lightblue", mid = "blue", high = "darkblue", midpoint = 0.7)
  
  # Save summary plot
  ggsave(
    filename = paste0(plotPath, "trajectory_correlations_", cellType, "_summary.png"),
    plot = summaryPlot,
    width = 10, height = max(6, nrow(retainedStats) * 0.4), dpi = 300
  )
  
  return(summaryPlot)
}
