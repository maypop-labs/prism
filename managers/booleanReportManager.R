# =============================================================================
# booleanReportManager.R (Cleaned Version - Only Functions Used by Script 06)
# Purpose: Boolean rule analysis and visualization for PRISM workflow
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)

# =============================================================================
# Main Report Generation Function
# =============================================================================

#' Generate comprehensive Boolean rule analysis report
#'
#' Creates enhanced visualizations and analysis incorporating SCENIC metadata
#' and multi-method rule inference results. Generates plots for rule quality
#' distribution, method comparison, and network complexity, plus text reports.
#'
#' @param boolRules List of Boolean rules from main script with rule metadata
#' @param edges Data frame with SCENIC metadata (TF, Target, corr, regType, etc.)
#' @param paths Paths object with directory structure containing base$plots and base$txt
#' @param cellType Character string with cell type name for file naming
#' @param trajectory Character string with trajectory name for file naming  
#' @param config Configuration object with plot parameters (figWidth, figHeight, figDPI)
#' @return Invisible list of generated plots, or NULL if no data available
#' @export
generateBooleanRuleReport <- function(boolRules, edges, paths, cellType, trajectory, config) {
  
  if (!dir.exists(paths$base$plots)) dir.create(paths$base$plots, recursive = TRUE)
  if (!dir.exists(paths$base$txt)) dir.create(paths$base$txt, recursive = TRUE)
  
  message("Generating SCENIC-enhanced Boolean rule report...")
  
  # Extract enhanced rule statistics
  ruleStats <- extractEnhancedRuleStatistics(boolRules, edges)
  
  if (nrow(ruleStats) == 0) {
    warning("No rule statistics available for reporting")
    return(invisible(NULL))
  }
  
  # Standard plot theme
  standardTheme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  # Generate enhanced visualizations
  plots <- list()
  
  tryCatch({
    # 1. Rule quality distribution with method breakdown
    plots$quality_dist <- plotRuleQualityByMethod(ruleStats, standardTheme)
    
    # 2. Method comparison
    plots$method_comparison <- plotMethodComparisonAnalysis(ruleStats, standardTheme)
    
    # 3. Network complexity overview
    plots$network_complexity <- plotNetworkComplexityOverview(ruleStats, standardTheme)
    
  }, error = function(e) {
    warning("Error generating plots: ", e$message)
  })
  
  # Save all plots
  plot_width <- config$figWidth %||% 10
  plot_height <- config$figHeight %||% 8
  plot_dpi <- config$figDPI %||% 300
  
  for (plot_name in names(plots)) {
    filename <- file.path(paths$base$plots, 
                         paste0(cellType, "_", trajectory, "_boolean_", plot_name, ".png"))
    
    tryCatch({
      ggsave(filename, plots[[plot_name]], 
             width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
      message("Saved plot: ", basename(filename))
    }, error = function(e) {
      warning("Failed to save plot ", plot_name, ": ", e$message)
    })
  }
  
  # Generate enhanced text reports
  tryCatch({
    generateEnhancedTextReport(ruleStats, boolRules, edges,
                              file.path(paths$base$txt, paste0(cellType, "_", trajectory, "_boolean_report.txt")))
  }, error = function(e) {
    warning("Error generating text reports: ", e$message)
  })
  
  message("Enhanced Boolean rule report generated successfully!")
  message("  - Plots saved to: ", paths$base$plots)
  message("  - Reports saved to: ", paths$base$txt)
  
  return(invisible(plots))
}

# =============================================================================
# Statistics Extraction
# =============================================================================

#' Extract comprehensive rule statistics with SCENIC metadata
#'
#' Processes Boolean rules to extract key statistics including scores, methods,
#' regulator counts, and rule complexity metrics. Categorizes rules by inference
#' method and analyzes logical operators used.
#'
#' @param boolRules List of Boolean rules with metadata (score, method, n_regulators, rule)
#' @param edges Data frame with SCENIC metadata (currently used for potential expansion)
#' @return Data frame with rule statistics and complexity metrics
#' @export
extractEnhancedRuleStatistics <- function(boolRules, edges) {
  
  if (length(boolRules) == 0) {
    return(data.frame())
  }
  
  # Basic rule statistics
  ruleStats <- data.frame(
    gene = names(boolRules),
    score = sapply(boolRules, function(x) x$score %||% 0),
    method = sapply(boolRules, function(x) x$method %||% "unknown"),
    n_regulators = sapply(boolRules, function(x) x$n_regulators %||% 0),
    rule_string = sapply(boolRules, function(x) x$rule %||% ""),
    stringsAsFactors = FALSE
  )
  
  # Add rule complexity metrics
  ruleStats$has_and <- grepl("&", ruleStats$rule_string)
  ruleStats$has_or <- grepl("\\|", ruleStats$rule_string)
  ruleStats$has_not <- grepl("!", ruleStats$rule_string)
  
  # Method categorization for analysis
  ruleStats$method_category <- sapply(ruleStats$method, function(m) {
    if (grepl("template", m)) return("Template-based")
    if (grepl("empirical", m)) return("Empirical")
    if (grepl("fallback|self_activation", m)) return("Fallback")
    return("Other")
  })
  
  return(ruleStats)
}

# =============================================================================
# Visualization Functions
# =============================================================================

#' Plot rule quality distribution by method
#'
#' Creates histogram showing distribution of rule quality scores faceted by
#' inference method category. Includes reference lines for quality thresholds.
#'
#' @param ruleStats Data frame with rule statistics from extractEnhancedRuleStatistics
#' @param standardTheme ggplot2 theme object for consistent plot styling
#' @return ggplot2 object with quality distribution histogram
#' @export
plotRuleQualityByMethod <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
           standardTheme +
           labs(title = "Rule Quality Distribution - No Data"))
  }
  
  ggplot(ruleStats, aes(x = score, fill = method_category)) +
    geom_histogram(bins = 25, alpha = 0.7, color = "white", size = 0.3) +
    facet_wrap(~method_category, scales = "free_y") +
    geom_vline(xintercept = 0.75, color = "red", linetype = "dashed", size = 0.5) +
    scale_fill_brewer(palette = "Set2", name = "Method") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Rule Quality Distribution by Method",
      subtitle = paste0("n = ", nrow(ruleStats), " rules across method categories"),
      x = "Rule Quality Score",
      y = "Number of Rules",
      caption = "Red line indicates high quality threshold (0.75)"
    ) +
    standardTheme +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"  # Remove legend since info is in facets
    )
}

#' Plot method comparison analysis
#'
#' Creates bar chart comparing mean performance across different rule inference
#' methods, showing sample sizes and relative performance.
#'
#' @param ruleStats Data frame with rule statistics from extractEnhancedRuleStatistics
#' @param standardTheme ggplot2 theme object for consistent plot styling
#' @return ggplot2 object with method comparison bar chart
#' @export
plotMethodComparisonAnalysis <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
           standardTheme +
           labs(title = "Method Comparison Analysis - No Data"))
  }
  
  # Method performance summary
  methodSummary <- ruleStats %>%
    group_by(method_category) %>%
    summarise(
      n_rules = n(),
      mean_score = mean(score, na.rm = TRUE),
      high_quality = sum(score >= 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      high_quality_pct = high_quality / n_rules * 100,
      method_category = reorder(method_category, mean_score)
    )
  
  ggplot(methodSummary, aes(x = method_category, y = mean_score)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0("n=", n_rules)), 
              vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Method Performance Comparison",
      subtitle = "Mean rule quality score by inference method",
      x = "Method Category",
      y = "Mean Rule Quality Score",
      caption = "Numbers above bars show sample size for each method"
    ) +
    standardTheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot network complexity overview
#'
#' Creates bar chart showing distribution of rules by number of regulators,
#' providing insight into network complexity and rule structure.
#'
#' @param ruleStats Data frame with rule statistics from extractEnhancedRuleStatistics
#' @param standardTheme ggplot2 theme object for consistent plot styling
#' @return ggplot2 object with complexity distribution bar chart
#' @export
plotNetworkComplexityOverview <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
           standardTheme +
           labs(title = "Network Complexity Overview - No Data"))
  }
  
  # Create complexity summaries
  regulator_dist <- data.frame(
    n_regulators = 0:max(ruleStats$n_regulators, na.rm = TRUE),
    count = sapply(0:max(ruleStats$n_regulators, na.rm = TRUE), 
                   function(x) sum(ruleStats$n_regulators == x, na.rm = TRUE))
  )
  
  ggplot(regulator_dist, aes(x = factor(n_regulators), y = count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 3) +
    labs(
      title = "Network Complexity: Rules by Number of Regulators",
      subtitle = paste0("Distribution across ", nrow(ruleStats), " total rules"),
      x = "Number of Regulators per Rule",
      y = "Number of Rules",
      caption = "Higher regulator counts indicate more complex regulatory logic"
    ) +
    standardTheme
}

# =============================================================================
# Text Reports
# =============================================================================

#' Generate enhanced text report with SCENIC integration
#'
#' Creates comprehensive text summary of Boolean rule inference results
#' including overall statistics, method performance breakdown, and key findings.
#'
#' @param ruleStats Data frame with rule statistics from extractEnhancedRuleStatistics
#' @param boolRules List of Boolean rules with metadata (used for additional analysis)
#' @param edges Data frame with SCENIC metadata (used for regulatory network context)
#' @param filename Character string with output file path for text report
#' @return NULL (writes file as side effect)
#' @export
generateEnhancedTextReport <- function(ruleStats, boolRules, edges, filename) {
  
  sink(filename)
  
  cat("=============================================================================\n")
  cat("BOOLEAN RULE INFERENCE REPORT (SCENIC-INTEGRATED)\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=============================================================================\n\n")
  
  # Overall Statistics
  cat("OVERALL STATISTICS:\n")
  cat("- Total genes with rules:", nrow(ruleStats), "\n")
  if (nrow(ruleStats) > 0) {
    cat("- Mean rule quality:", round(mean(ruleStats$score, na.rm = TRUE), 3), "\n")
    cat("- Median rule quality:", round(median(ruleStats$score, na.rm = TRUE), 3), "\n")
    cat("- Rules with quality >= 0.75:", sum(ruleStats$score >= 0.75, na.rm = TRUE), 
        "(", round(100 * sum(ruleStats$score >= 0.75, na.rm = TRUE) / nrow(ruleStats), 1), "%)")
    cat("- Rules with quality >= 0.9:", sum(ruleStats$score >= 0.9, na.rm = TRUE), 
        "(", round(100 * sum(ruleStats$score >= 0.9, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n\n")
  }
  
  # Method Analysis
  if (nrow(ruleStats) > 0) {
    cat("METHOD PERFORMANCE:\n")
    methodCounts <- table(ruleStats$method_category)
    for (method in names(methodCounts)) {
      method_rules <- ruleStats[ruleStats$method_category == method, ]
      cat("- ", method, ":", methodCounts[method], "rules\n")
      cat("  • Mean quality:", round(mean(method_rules$score, na.rm = TRUE), 3), "\n")
      cat("  • High quality (>=0.75):", sum(method_rules$score >= 0.75, na.rm = TRUE), 
          "/", nrow(method_rules), "\n")
    }
    cat("\n")
  }
  
  # Network Complexity
  if (nrow(ruleStats) > 0) {
    cat("NETWORK COMPLEXITY:\n")
    cat("- Genes with 0 regulators:", sum(ruleStats$n_regulators == 0, na.rm = TRUE), "\n")
    cat("- Genes with 1 regulator:", sum(ruleStats$n_regulators == 1, na.rm = TRUE), "\n")
    cat("- Genes with 2 regulators:", sum(ruleStats$n_regulators == 2, na.rm = TRUE), "\n")
    cat("- Genes with 3+ regulators:", sum(ruleStats$n_regulators >= 3, na.rm = TRUE), "\n")
    cat("- Mean regulators per gene:", round(mean(ruleStats$n_regulators, na.rm = TRUE), 2), "\n")
    cat("- Rules with AND logic:", sum(ruleStats$has_and), "\n")
    cat("- Rules with OR logic:", sum(ruleStats$has_or), "\n")
    cat("- Rules with NOT logic:", sum(ruleStats$has_not), "\n\n")
  }
  
  # Top Quality Rules
  if (nrow(ruleStats) > 0) {
    cat("TOP 10 HIGHEST QUALITY RULES:\n")
    topRules <- ruleStats[order(ruleStats$score, decreasing = TRUE, na.last = TRUE)[1:min(10, nrow(ruleStats))], ]
    for (i in 1:nrow(topRules)) {
      cat(i, ".", topRules$gene[i], "- Score:", round(topRules$score[i], 3), 
          "- Method:", topRules$method[i], "- Regulators:", topRules$n_regulators[i], "\n")
    }
    cat("\n")
  }
  
  # SCENIC Edge Analysis
  if (nrow(edges) > 0) {
    cat("SCENIC REGULATORY NETWORK CONTEXT:\n")
    cat("- Total regulatory edges:", nrow(edges), "\n")
    cat("- Unique transcription factors:", length(unique(edges$TF)), "\n")
    cat("- Unique target genes:", length(unique(edges$Target)), "\n")
    
    if ("regType" %in% colnames(edges)) {
      regTypeCounts <- table(edges$regType, useNA = "ifany")
      for (type in names(regTypeCounts)) {
        type_label <- if (is.na(type)) "Unknown" else type
        cat("- ", type_label, "edges:", regTypeCounts[type], "\n")
      }
    }
    
    if ("hasMotif" %in% colnames(edges)) {
      cat("- Edges with motif support:", sum(edges$hasMotif %in% TRUE, na.rm = TRUE), 
          "/", nrow(edges), "\n")
    }
    cat("\n")
  }
  
  cat("=============================================================================\n")
  cat("END OF REPORT\n")
  cat("=============================================================================\n")
  
  sink()
  
  message("Enhanced text report saved to: ", filename)
}

# =============================================================================
# Utility Functions
# =============================================================================

#' Null coalescing operator
#'
#' Returns the right-hand side if left-hand side is NULL, empty, or NA.
#'
#' @param x Left-hand side value
#' @param y Right-hand side fallback value
#' @return x if valid, otherwise y
#' @export
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x
