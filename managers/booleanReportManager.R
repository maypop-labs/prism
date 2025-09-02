# =============================================================================
# booleanReportManager.R (Production Version with SCENIC Integration)
# Purpose: SCENIC-enhanced Boolean rule analysis and visualization
# Enhanced with: regulator weight analysis, sign consistency metrics, method comparison
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
#' and multi-method rule inference results.
#'
#' @param boolRules List of Boolean rules from main script
#' @param edges Edge list data frame with SCENIC metadata
#' @param paths Paths object with directory structure
#' @param cellType Cell type name for file naming
#' @param trajectory Trajectory name for file naming  
#' @param config Configuration object with plot parameters
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
    
    # 2. SCENIC integration analysis
    plots$scenic_analysis <- plotScenicIntegrationAnalysis(ruleStats, standardTheme)
    
    # 3. Sign consistency analysis
    plots$sign_consistency <- plotSignConsistencyAnalysis(ruleStats, standardTheme)
    
    # 4. Regulator weight distribution
    plots$regulator_weights <- plotRegulatorWeightDistribution(ruleStats, edges, standardTheme)
    
    # 5. Method comparison
    plots$method_comparison <- plotMethodComparisonAnalysis(ruleStats, standardTheme)
    
    # 6. Network complexity overview
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
    
    generateSummaryStatistics(ruleStats, boolRules, edges,
                             file.path(paths$base$txt, paste0(cellType, "_", trajectory, "_boolean_summary.txt")))
  }, error = function(e) {
    warning("Error generating text reports: ", e$message)
  })
  
  message("Enhanced Boolean rule report generated successfully!")
  message("  - Plots saved to: ", paths$base$plots)
  message("  - Reports saved to: ", paths$base$txt)
  
  return(invisible(plots))
}

# =============================================================================
# Enhanced Statistics Extraction
# =============================================================================

#' Extract comprehensive rule statistics with SCENIC metadata
#'
#' @param boolRules List of Boolean rules
#' @param edges Edge list with SCENIC metadata
#' @return Data frame with enhanced statistics
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
    k_used = sapply(boolRules, function(x) x$k_used %||% 0),
    stringsAsFactors = FALSE
  )
  
  # Add rule complexity metrics
  ruleStats$has_and <- grepl("&", ruleStats$rule_string)
  ruleStats$has_or <- grepl("\\|", ruleStats$rule_string)
  ruleStats$has_not <- grepl("!", ruleStats$rule_string)
  ruleStats$rule_length <- nchar(ruleStats$rule_string)
  
  # Count operators more carefully
  ruleStats$n_operators <- sapply(ruleStats$rule_string, function(rule) {
    logic_part <- strsplit(rule, ",\\s*")[[1]]
    if (length(logic_part) >= 2) {
      logic <- logic_part[2]
      and_count <- length(gregexpr("&", logic, fixed = TRUE)[[1]])
      or_count <- length(gregexpr("\\|", logic)[[1]])
      not_count <- length(gregexpr("!", logic, fixed = TRUE)[[1]])
      
      # Adjust for no matches (-1 return)
      and_count <- ifelse(and_count == 1 && !grepl("&", logic), 0, and_count)
      or_count <- ifelse(or_count == 1 && !grepl("\\|", logic), 0, or_count)
      not_count <- ifelse(not_count == 1 && !grepl("!", logic), 0, not_count)
      
      return(and_count + or_count + not_count)
    }
    return(0)
  })
  
  # Add score breakdown if available
  ruleStats$base_accuracy <- sapply(boolRules, function(x) {
    if (!is.null(x$score_breakdown)) x$score_breakdown$base_accuracy else x$score %||% 0
  })
  
  ruleStats$sign_consistency <- sapply(boolRules, function(x) {
    if (!is.null(x$score_breakdown)) x$score_breakdown$sign_consistency else 0
  })
  
  ruleStats$scenic_confidence <- sapply(boolRules, function(x) {
    if (!is.null(x$score_breakdown)) x$score_breakdown$scenic_confidence else 0
  })
  
  # Method categorization
  ruleStats$method_category <- sapply(ruleStats$method, function(m) {
    if (grepl("template", m)) return("Template-based")
    if (grepl("empirical", m)) return("Empirical")
    if (grepl("fallback|self_activation", m)) return("Fallback")
    return("Other")
  })
  
  # Add regulator information from edges
  if (nrow(edges) > 0) {
    ruleStats$has_scenic_regs <- sapply(ruleStats$gene, function(g) {
      target_edges <- edges[edges$Target == g, , drop = FALSE]
      nrow(target_edges) > 0
    })
    
    ruleStats$mean_reg_corr <- sapply(ruleStats$gene, function(g) {
      target_edges <- edges[edges$Target == g, , drop = FALSE]
      if (nrow(target_edges) > 0) mean(abs(target_edges$corr), na.rm = TRUE) else 0
    })
  } else {
    ruleStats$has_scenic_regs <- FALSE
    ruleStats$mean_reg_corr <- 0
  }
  
  return(ruleStats)
}

# =============================================================================
# Enhanced Visualization Functions
# =============================================================================

#' Plot rule quality distribution by method
plotRuleQualityByMethod <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
           standardTheme +
           labs(title = "Rule Quality Distribution - No Data"))
  }
  
  # Create quality categories
  ruleStats$quality_category <- cut(ruleStats$score, 
                                    breaks = c(-Inf, 0.5, 0.75, 0.9, Inf),
                                    labels = c("Low (<0.5)", "Medium (0.5-0.75)", 
                                               "High (0.75-0.9)", "Excellent (>0.9)"))
  
  ggplot(ruleStats, aes(x = score, fill = method_category)) +
    geom_histogram(bins = 25, alpha = 0.7, color = "white", size = 0.3) +
    facet_wrap(~method_category, scales = "free_y") +
    geom_vline(xintercept = 0.75, color = "red", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 0.9, color = "darkgreen", linetype = "dashed", size = 0.5) +
    scale_fill_brewer(palette = "Set2", name = "Method") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Rule Quality Distribution by Method",
      subtitle = paste0("n = ", nrow(ruleStats), " rules across ", 
                        length(unique(ruleStats$method_category)), " method categories"),
      x = "Rule Quality Score",
      y = "Number of Rules",
      caption = "Red line: good threshold (0.75), Green line: excellent threshold (0.9)"
    ) +
    standardTheme +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"  # Remove legend since info is in facets
    )
}

#' Plot SCENIC integration analysis
plotScenicIntegrationAnalysis <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0 || !("scenic_confidence" %in% colnames(ruleStats))) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No SCENIC data available", size = 5) +
           standardTheme +
           labs(title = "SCENIC Integration Analysis - No Data"))
  }
  
  # Only analyze rules with SCENIC regulators
  scenic_rules <- ruleStats[ruleStats$has_scenic_regs, ]
  
  if (nrow(scenic_rules) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No rules with SCENIC regulators", size = 5) +
           standardTheme +
           labs(title = "SCENIC Integration Analysis - No SCENIC Rules"))
  }
  
  # Calculate correlation
  correlation <- cor(scenic_rules$scenic_confidence, scenic_rules$score, 
                     use = "complete.obs", method = "pearson")
  
  ggplot(scenic_rules, aes(x = scenic_confidence, y = score, color = method_category)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3, color = "black") +
    scale_color_brewer(palette = "Set1", name = "Method") +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Rule Quality vs SCENIC Confidence",
      subtitle = "Higher SCENIC confidence should correlate with better rule quality",
      x = "SCENIC Confidence Score",
      y = "Rule Quality Score",
      caption = paste0("Correlation: r = ", round(correlation, 3), " (n = ", nrow(scenic_rules), " rules)")
    ) +
    standardTheme
}

#' Plot sign consistency analysis
plotSignConsistencyAnalysis <- function(ruleStats, standardTheme) {
  
  if (nrow(ruleStats) == 0 || !("sign_consistency" %in% colnames(ruleStats))) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No sign consistency data available", size = 5) +
           standardTheme +
           labs(title = "Sign Consistency Analysis - No Data"))
  }
  
  # Filter to rules that have regulators
  regRules <- ruleStats[ruleStats$n_regulators > 0, ]
  
  if (nrow(regRules) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No rules with regulators", size = 5) +
           standardTheme +
           labs(title = "Sign Consistency Analysis - No Regulator Rules"))
  }
  
  ggplot(regRules, aes(x = sign_consistency, fill = method_category)) +
    geom_histogram(bins = 20, alpha = 0.7, color = "white") +
    scale_fill_brewer(palette = "Set2", name = "Method") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Sign Consistency Analysis",
      subtitle = paste0("Distribution of sign consistency scores (n = ", nrow(regRules), " rules with regulators)"),
      x = "Sign Consistency Score (0 = poor, 1 = perfect)",
      y = "Number of Rules",
      caption = "Measures how well inferred rules match SCENIC regulatory directions"
    ) +
    standardTheme +
    facet_wrap(~method_category, scales = "free_y") +
    theme(strip.text = element_text(face = "bold"))
}

#' Plot regulator weight distribution
plotRegulatorWeightDistribution <- function(ruleStats, edges, standardTheme) {
  
  if (nrow(edges) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "No edge data available", size = 5) +
           standardTheme +
           labs(title = "Regulator Weight Distribution - No Data"))
  }
  
  # Calculate weights for all edges (using same logic as booleanManager)
  edges$abs_corr <- abs(edges$corr)
  edges$motif_bonus <- ifelse(!is.na(edges$motifConfidence), edges$motifConfidence, 0)
  edges$nes_bonus <- ifelse(!is.na(edges$NES), pmax(0, edges$NES / 5), 0)
  edges$prior_bonus <- ifelse(!is.na(edges$hasMotif) & edges$hasMotif == TRUE, 0.1, 0)
  
  edges$composite_weight <- edges$abs_corr + 
                           (edges$motif_bonus * 0.3) + 
                           (edges$nes_bonus * 0.2) + 
                           edges$prior_bonus
  
  # Separate used vs unused regulators
  used_targets <- ruleStats$gene[ruleStats$n_regulators > 0]
  edges$used_in_rules <- edges$Target %in% used_targets
  
  ggplot(edges, aes(x = composite_weight, fill = used_in_rules)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "steelblue"),
                      name = "Used in Rules",
                      labels = c("Not Used", "Used")) +
    labs(
      title = "Regulator Weight Distribution",
      subtitle = paste0("Composite weights for ", nrow(edges), " regulatory relationships"),
      x = "Composite Weight (correlation + motif + NES + prior)",
      y = "Number of Regulatory Edges",
      caption = "Higher weights indicate stronger evidence for regulatory relationship"
    ) +
    standardTheme
}

#' Plot method comparison analysis
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
      median_score = median(score, na.rm = TRUE),
      high_quality = sum(score >= 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      high_quality_pct = high_quality / n_rules * 100,
      method_category = reorder(method_category, mean_score)
    )
  
  p1 <- ggplot(methodSummary, aes(x = method_category, y = mean_score)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0("n=", n_rules)), 
              vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Method Performance Comparison",
      x = "Method Category",
      y = "Mean Rule Quality Score"
    ) +
    standardTheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(methodSummary, aes(x = method_category, y = high_quality_pct)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    geom_text(aes(label = paste0(round(high_quality_pct, 1), "%")), 
              vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      title = "High Quality Rules by Method",
      x = "Method Category", 
      y = "% Rules with Score ≥ 0.75"
    ) +
    standardTheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2, top = "Method Comparison Analysis")
}

#' Plot network complexity overview
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
  
  logic_summary <- data.frame(
    logic_type = c("AND only", "OR only", "NOT only", "Mixed", "Simple"),
    count = c(
      sum(ruleStats$has_and & !ruleStats$has_or & !ruleStats$has_not),
      sum(!ruleStats$has_and & ruleStats$has_or & !ruleStats$has_not),
      sum(!ruleStats$has_and & !ruleStats$has_or & ruleStats$has_not),
      sum((ruleStats$has_and + ruleStats$has_or + ruleStats$has_not) >= 2),
      sum(!ruleStats$has_and & !ruleStats$has_or & !ruleStats$has_not)
    )
  )
  
  p1 <- ggplot(regulator_dist, aes(x = factor(n_regulators), y = count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 3) +
    labs(
      title = "Rules by Number of Regulators",
      x = "Number of Regulators",
      y = "Number of Rules"
    ) +
    standardTheme
  
  p2 <- ggplot(logic_summary[logic_summary$count > 0, ], aes(x = reorder(logic_type, count), y = count)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    geom_text(aes(label = count), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(
      title = "Rules by Logic Type",
      x = "Logic Type",
      y = "Number of Rules"
    ) +
    standardTheme
  
  gridExtra::grid.arrange(p1, p2, ncol = 2, top = "Network Complexity Overview")
}

# =============================================================================
# Enhanced Text Reports
# =============================================================================

#' Generate enhanced text report with SCENIC integration
generateEnhancedTextReport <- function(ruleStats, boolRules, edges, filename) {
  
  sink(filename)
  
  cat("=============================================================================\n")
  cat("ENHANCED BOOLEAN RULE INFERENCE REPORT (SCENIC-INTEGRATED)\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=============================================================================\n\n")
  
  # Overall Statistics
  cat("OVERALL STATISTICS:\n")
  cat("- Total genes with rules:", nrow(ruleStats), "\n")
  if (nrow(ruleStats) > 0) {
    cat("- Mean rule quality:", round(mean(ruleStats$score, na.rm = TRUE), 3), "\n")
    cat("- Median rule quality:", round(median(ruleStats$score, na.rm = TRUE), 3), "\n")
    cat("- Rules with quality ≥ 0.75:", sum(ruleStats$score >= 0.75, na.rm = TRUE), 
        "(", round(100 * sum(ruleStats$score >= 0.75, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n")
    cat("- Rules with quality ≥ 0.9:", sum(ruleStats$score >= 0.9, na.rm = TRUE), 
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
      cat("  • High quality (≥0.75):", sum(method_rules$score >= 0.75, na.rm = TRUE), 
          "/", nrow(method_rules), "\n")
    }
    cat("\n")
  }
  
  # SCENIC Integration Analysis
  if (nrow(ruleStats) > 0 && "scenic_confidence" %in% colnames(ruleStats)) {
    cat("SCENIC INTEGRATION ANALYSIS:\n")
    cat("- Rules with SCENIC regulators:", sum(ruleStats$has_scenic_regs), "/", nrow(ruleStats), "\n")
    cat("- Mean SCENIC confidence:", round(mean(ruleStats$scenic_confidence, na.rm = TRUE), 3), "\n")
    
    if ("sign_consistency" %in% colnames(ruleStats)) {
      reg_rules <- ruleStats[ruleStats$n_regulators > 0, ]
      if (nrow(reg_rules) > 0) {
        cat("- Mean sign consistency:", round(mean(reg_rules$sign_consistency, na.rm = TRUE), 3), "\n")
        cat("- Perfect sign consistency:", sum(reg_rules$sign_consistency >= 0.99, na.rm = TRUE), 
            "/", nrow(reg_rules), "\n")
      }
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
    cat("SCENIC EDGE ANALYSIS:\n")
    cat("- Total regulatory edges:", nrow(edges), "\n")
    cat("- Unique TFs:", length(unique(edges$TF)), "\n")
    cat("- Unique targets:", length(unique(edges$Target)), "\n")
    
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
  
  sink()
  
  message("Enhanced text report saved to: ", filename)
}

#' Generate summary statistics file
generateSummaryStatistics <- function(ruleStats, boolRules, edges, filename) {
  
  sink(filename)
  
  cat("BOOLEAN RULE INFERENCE SUMMARY STATISTICS\n")
  cat("=========================================\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  if (nrow(ruleStats) > 0) {
    cat("RULE QUALITY METRICS:\n")
    cat("Count:", nrow(ruleStats), "\n")
    cat("Mean:", round(mean(ruleStats$score, na.rm = TRUE), 4), "\n")
    cat("Median:", round(median(ruleStats$score, na.rm = TRUE), 4), "\n")
    cat("SD:", round(sd(ruleStats$score, na.rm = TRUE), 4), "\n")
    cat("Min:", round(min(ruleStats$score, na.rm = TRUE), 4), "\n")
    cat("Max:", round(max(ruleStats$score, na.rm = TRUE), 4), "\n")
    cat("Q1:", round(quantile(ruleStats$score, 0.25, na.rm = TRUE), 4), "\n")
    cat("Q3:", round(quantile(ruleStats$score, 0.75, na.rm = TRUE), 4), "\n\n")
    
    cat("QUALITY THRESHOLDS:\n")
    cat("≥0.5:", sum(ruleStats$score >= 0.5, na.rm = TRUE), "/", nrow(ruleStats), 
        "(", round(100 * mean(ruleStats$score >= 0.5, na.rm = TRUE), 1), "%)\n")
    cat("≥0.75:", sum(ruleStats$score >= 0.75, na.rm = TRUE), "/", nrow(ruleStats), 
        "(", round(100 * mean(ruleStats$score >= 0.75, na.rm = TRUE), 1), "%)\n")
    cat("≥0.9:", sum(ruleStats$score >= 0.9, na.rm = TRUE), "/", nrow(ruleStats), 
        "(", round(100 * mean(ruleStats$score >= 0.9, na.rm = TRUE), 1), "%)\n\n")
  }
  
  if (nrow(edges) > 0) {
    cat("REGULATORY NETWORK METRICS:\n")
    cat("Total edges:", nrow(edges), "\n")
    cat("Unique TFs:", length(unique(edges$TF)), "\n")
    cat("Unique targets:", length(unique(edges$Target)), "\n")
    cat("Mean |correlation|:", round(mean(abs(edges$corr), na.rm = TRUE), 4), "\n")
    
    if ("motifConfidence" %in% colnames(edges)) {
      cat("Mean motif confidence:", round(mean(edges$motifConfidence, na.rm = TRUE), 4), "\n")
    }
    if ("NES" %in% colnames(edges)) {
      cat("Mean NES score:", round(mean(edges$NES, na.rm = TRUE), 4), "\n")
    }
  }
  
  sink()
  
  message("Summary statistics saved to: ", filename)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x
