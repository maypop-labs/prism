# =============================================================================
# booleanReportManager.R
# Purpose: Boolean network plotting and report file generation.
# 
# =============================================================================

# =============================================================================
# Boolean Rule Analysis and Visualization Functions
# =============================================================================

#' Generate comprehensive Boolean rule analysis report
#'
#' Creates visualizations using config parameters and proper white backgrounds
#' Called by: Main script(s)!
#'
#' @param boolRules List of Boolean rules from main script
#' @param edges Edge list data frame with regulatory relationships
#' @param paths Paths object with directory structure
#' @param cellType Cell type name for file naming
#' @param trajectory Trajectory name for file naming
#' @param config Configuration object with plot parameters
#' @export
generateBooleanRuleReport <- function(boolRules, edges, paths, cellType, trajectory, config) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(ComplexHeatmap)
  library(igraph)
  
  if (!dir.exists(paths$base$plots)) dir.create(paths$base$plots, recursive = TRUE)
  if (!dir.exists(paths$base$txt)) dir.create(paths$base$txt, recursive = TRUE)
  
  # Extract summary statistics
  ruleStats <- extractRuleStatistics(boolRules)
  
  # Standard plot theme using config parameters
  standardTheme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  # 1. Rule Quality Distribution
  p1 <- plotRuleQualityDistribution(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_rule_quality_distribution.png"), 
         p1, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 2. Method Usage Summary
  p2 <- plotMethodUsageSummary(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_method_usage_summary.png"), 
         p2, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 3. Biological Plausibility Analysis
  p3 <- plotBiologicalPlausibility(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_biological_plausibility.png"), 
         p3, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 4. Network Complexity Heatmap (FIXED VERSION)
  if (length(boolRules) > 0) {
    p4 <- plotNetworkComplexityHeatmap(boolRules, edges)
    png(paste0(paths$base$plots, cellType, "_", trajectory, "_network_complexity_heatmap.png"), 
        width = config$figWidth * config$figDPI, height = config$figHeight * config$figDPI, 
        res = config$figDPI, bg = "white")
    draw(p4)
    dev.off()
  }
  
  # 5. Rule Logic Summary Network
  p5 <- plotRuleLogicNetwork(boolRules, edges, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_rule_logic_network.png"), 
         p5, width = config$figWidth * 1.5, height = config$figHeight * 1.2, dpi = config$figDPI, bg = "white")
  
  # 6. Generate text report
  generateTextReport(ruleStats, boolRules, paste0(paths$base$txt, cellType, "_", trajectory, "_boolean_rules_report.txt"))
  
  message("Boolean rule analysis complete. Results saved to: ", paths$base$plots)
}

#' Plot distribution of rule quality scores
#' Called by: generateBooleanRuleReport
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotRuleQualityDistribution <- function(ruleStats, standardTheme) {
  
  ggplot(ruleStats, aes(x = bestScore)) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.95, color = "darkgreen", linetype = "dashed", size = 1) +
    labs(
      title = "Boolean Rule Quality Distribution",
      subtitle = paste0("n = ", nrow(ruleStats), " genes with inferred rules"),
      x = "Rule Quality Score (Accuracy)",
      y = "Number of Genes",
      caption = "Red line: minimum threshold (0.8), Green line: excellent threshold (0.95)"
    ) +
    standardTheme +
    annotate("text", x = 0.85, y = Inf, label = "Good Rules", vjust = 1.5, color = "darkgreen", size = 3) +
    annotate("text", x = 0.4, y = Inf, label = "Poor Rules", vjust = 1.5, color = "red", size = 3)
}

#' Plot method usage summary
#' Called by: generateBooleanRuleReport
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotMethodUsageSummary <- function(ruleStats, standardTheme) {
  
  methodCounts <- ruleStats %>%
    dplyr::group_by(methodUsed) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::mutate(
      percentage = round(100 * count / sum(count), 1),
      label = paste0(percentage, "%")
    )
  
  ggplot(methodCounts, aes(x = "", y = count, fill = methodUsed)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(
      title = "Boolean Rule Inference Method Usage",
      subtitle = "Which methods were used to find the best rules?",
      fill = "Method"
    ) +
    standardTheme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
}

#' Plot biological plausibility analysis
#' Called by: generateBooleanRuleReport
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotBiologicalPlausibility <- function(ruleStats, standardTheme) {
  
  plausibilityData <- ruleStats %>%
    filter(!is.na(biologicallyPlausible)) %>%
    mutate(
      plausibilityLabel = ifelse(biologicallyPlausible, "Biologically Plausible", "Data-Driven Only")
    )
  
  if (nrow(plausibilityData) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No biological plausibility data available", size = 5) +
             labs(title = "Biological Plausibility Analysis", subtitle = "No data available") +
             standardTheme)
  }
  
  ggplot(plausibilityData, aes(x = plausibilityLabel, y = bestScore, fill = plausibilityLabel)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    labs(
      title = "Rule Quality by Biological Plausibility",
      subtitle = "Do biologically plausible rules perform better?",
      x = "Rule Type",
      y = "Quality Score",
      fill = "Rule Type"
    ) +
    standardTheme +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("steelblue", "orange")) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, color = "red", fill = "white")
}

#' Create network complexity heatmap
#' Called by: generateBooleanRuleReport
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @return ComplexHeatmap object
plotNetworkComplexityHeatmap <- function(boolRules, edges) {
  
  # FIXED: Only use genes that actually have Boolean rules and their regulators
  targetGenes <- names(boolRules)
  regulatorGenes <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
  relevantGenes <- unique(c(targetGenes, regulatorGenes))
  
  # Create appropriately sized matrix
  regulatorMatrix <- matrix(0, nrow = length(regulatorGenes), ncol = length(targetGenes))
  rownames(regulatorMatrix) <- regulatorGenes
  colnames(regulatorMatrix) <- targetGenes
  
  # Fill in regulatory relationships with quality scores
  for (gene in targetGenes) {
    if (gene %in% names(boolRules)) {
      regulators <- boolRules[[gene]]$regulators
      
      # Ensure ruleQuality is a single numeric value
      ruleQuality <- as.numeric(boolRules[[gene]]$bestScore %||% boolRules[[gene]]$score %||% 0)[1]
      if (is.na(ruleQuality)) ruleQuality <- 0
      
      for (reg in regulators) {
        if (reg %in% rownames(regulatorMatrix)) {
          regulatorMatrix[reg, gene] <- ruleQuality
        }
      }
    }
  }
  
  # Remove empty rows and columns
  nonEmptyRows <- rowSums(regulatorMatrix) > 0
  nonEmptyCols <- colSums(regulatorMatrix) > 0
  
  if (sum(nonEmptyRows) > 0 && sum(nonEmptyCols) > 0) {
    regulatorMatrix <- regulatorMatrix[nonEmptyRows, nonEmptyCols, drop = FALSE]
    
    # Create properly sized heatmap
    ComplexHeatmap::Heatmap(
      regulatorMatrix,
      name = "Rule Quality",
      col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
      cluster_rows = ifelse(nrow(regulatorMatrix) > 1, TRUE, FALSE),
      cluster_columns = ifelse(ncol(regulatorMatrix) > 1, TRUE, FALSE),
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10),
      column_title = "Target Genes",
      row_title = "Regulator Genes",
      heatmap_legend_param = list(title = "Boolean Rule\nQuality Score"),
      width = unit(min(8, ncol(regulatorMatrix) * 0.8), "cm"),
      height = unit(min(8, nrow(regulatorMatrix) * 0.8), "cm")
    )
  } else {
    # Create empty heatmap with message
    emptyMatrix <- matrix(0, nrow = 1, ncol = 1)
    rownames(emptyMatrix) <- "No Data"
    colnames(emptyMatrix) <- "No Data"
    
    ComplexHeatmap::Heatmap(
      emptyMatrix,
      name = "Rule Quality",
      col = "white",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = "No regulatory relationships to display",
      heatmap_legend_param = list(title = "Boolean Rule\nQuality Score")
    )
  }
}

#' Plot rule logic network
#' Called by: generateBooleanRuleReport
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotRuleLogicNetwork <- function(boolRules, edges, standardTheme) {
  
  # Filter edges to only include those with Boolean rules
  networkEdges <- edges %>%
    dplyr::filter(Target %in% names(boolRules)) %>%
    dplyr::mutate(
      ruleQuality = sapply(Target, function(x) {
        if (x %in% names(boolRules)) boolRules[[x]]$bestScore else 0
      })
    )
  
  if (nrow(networkEdges) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No network data to display", size = 5) +
             labs(title = "Gene Regulatory Network", subtitle = "No edges available") +
             standardTheme)
  }
  
  # Create igraph object
  g <- igraph::graph_from_data_frame(networkEdges, directed = TRUE)
  
  # Calculate layout
  if (igraph::vcount(g) > 1) {
    layout <- igraph::layout_with_fr(g)
  } else {
    layout <- matrix(c(0, 0), nrow = 1)
  }
  
  # Convert to data frame for ggplot
  nodeList <- data.frame(
    name = igraph::V(g)$name,
    x = layout[,1],
    y = layout[,2],
    isTarget = igraph::V(g)$name %in% names(boolRules)
  )
  
  # Create edge coordinates
  edgeList <- as.data.frame(igraph::get.edgelist(g))
  names(edgeList) <- c("from", "to")
  edgeList$quality <- igraph::E(g)$ruleQuality
  
  edgeCoords <- edgeList %>%
    dplyr::left_join(nodeList, by = c("from" = "name")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::select(-isTarget) %>%
    dplyr::left_join(nodeList, by = c("to" = "name")) %>%
    dplyr::rename(x2 = x, y2 = y, isTarget = isTarget)
  
  # Plot with better styling
  ggplot() +
    geom_segment(
      data = edgeCoords,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = quality),
      arrow = arrow(length = unit(0.15, "inches")),
      alpha = 0.7,
      size = 1
    ) +
    geom_point(
      data = nodeList,
      aes(x = x, y = y, shape = isTarget, fill = isTarget),
      size = 4,
      alpha = 0.8,
      color = "black",
      stroke = 0.5
    ) +
    scale_color_gradient(
      low = "lightgray",
      high = "red",
      name = "Rule\nQuality"
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 22),
      name = "Node Type",
      labels = c("Regulator", "Target")
    ) +
    scale_fill_manual(
      values = c("FALSE" = "lightblue", "TRUE" = "orange"),
      name = "Node Type",
      labels = c("Regulator", "Target")
    ) +
    labs(
      title = "Gene Regulatory Network with Boolean Rules",
      subtitle = paste0("Network includes ", nrow(nodeList), " genes and ", nrow(edgeCoords), " regulatory relationships")
    ) +
    standardTheme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

#' Extract summary statistics from Boolean rule list
#' Called by: generateBooleanRuleReport
#' @param boolRules List of Boolean rule objects
#' @return Data frame with per-gene statistics
extractRuleStatistics <- function(boolRules) {
  data.frame(
    gene = names(boolRules),
    nRegulators = sapply(boolRules, function(x) x$nRegulators %||% NA),
    bestScore = sapply(boolRules, function(x) x$score %||% NA),
    methodUsed = sapply(boolRules, function(x) x$methodUsed %||% "template"),
    biologicallyPlausible = sapply(boolRules, function(x) x$biologicallyPlausible %||% NA)
  )
}

#' Generate text summary report
#' Called by: generateBooleanRuleReport
#' @param ruleStats Data frame with rule statistics
#' @param boolRules List of Boolean rules
#' @param filename Output file path
generateTextReport <- function(ruleStats, boolRules, filename) {
  
  sink(filename)
  
  cat("=============================================================================\n")
  cat("BOOLEAN RULE INFERENCE SUMMARY REPORT\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=============================================================================\n\n")
  
  cat("OVERALL STATISTICS:\n")
  cat("- Total genes with rules:", nrow(ruleStats), "\n")
  cat("- Mean rule quality:", round(mean(ruleStats$bestScore, na.rm = TRUE), 3), "\n")
  cat("- Median rule quality:", round(median(ruleStats$bestScore, na.rm = TRUE), 3), "\n")
  cat("- Rules with quality > 0.8:", sum(ruleStats$bestScore > 0.8, na.rm = TRUE), 
      "(", round(100 * sum(ruleStats$bestScore > 0.8, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n")
  cat("- Rules with quality > 0.95:", sum(ruleStats$bestScore > 0.95, na.rm = TRUE), 
      "(", round(100 * sum(ruleStats$bestScore > 0.95, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n\n")
  
  cat("NETWORK COMPLEXITY:\n")
  cat("- Genes with 1 regulator:", sum(ruleStats$nRegulators == 1, na.rm = TRUE), "\n")
  cat("- Genes with 2 regulators:", sum(ruleStats$nRegulators == 2, na.rm = TRUE), "\n")
  cat("- Genes with 3+ regulators:", sum(ruleStats$nRegulators >= 3, na.rm = TRUE), "\n")
  cat("- Mean regulators per gene:", round(mean(ruleStats$nRegulators, na.rm = TRUE), 2), "\n\n")
  
  if ("methodUsed" %in% colnames(ruleStats)) {
    cat("METHOD USAGE:\n")
    methodCounts <- table(ruleStats$methodUsed, useNA = "no")
    for (method in names(methodCounts)) {
      cat("-", method, ":", methodCounts[method], "genes\n")
    }
    cat("\n")
  }
  
  cat("TOP 10 HIGHEST QUALITY RULES:\n")
  topRules <- ruleStats[order(ruleStats$bestScore, decreasing = TRUE, na.last = TRUE)[1:min(10, nrow(ruleStats))], ]
  for (i in 1:nrow(topRules)) {
    cat(i, ".", topRules$gene[i], "- Score:", round(topRules$bestScore[i], 3), 
        "- Regulators:", topRules$nRegulators[i], "\n")
  }
  
  sink()
  
  message("Text report saved to: ", filename)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x
