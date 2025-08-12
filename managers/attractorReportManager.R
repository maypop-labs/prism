# =============================================================================
# BoolNet Analysis Validation Functions
# =============================================================================

#' Quick validation summary of BoolNet results
#'
#' @param attractors Attractors object from getAttractors()
#' @param boolnet BoolNet network object
#' @return Summary statistics
validateBoolNetResults <- function(attractors, boolnet) {
  
  cat("=== BOOLNET ANALYSIS VALIDATION ===\n")
  
  # Basic network info
  cat("NETWORK:\n")
  cat("  Genes in network:", length(boolnet$genes), "\n")
  cat("  Network type:", boolnet$type, "\n\n")
  
  # Attractor summary
  nAttractors <- length(attractors$attractors)
  cat("ATTRACTORS:\n")
  cat("  Total attractors found:", nAttractors, "\n")
  
  if (nAttractors > 0) {
    # Attractor sizes
    attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))
    cat("  Attractor size range:", min(attractorSizes), "to", max(attractorSizes), "states\n")
    cat("  Steady states (size=1):", sum(attractorSizes == 1), "\n")
    cat("  Cycles (size>1):", sum(attractorSizes > 1), "\n")
    
    # Basin sizes
    if ("basinSize" %in% names(attractors$attractors[[1]])) {
      basinSizes <- sapply(attractors$attractors, function(att) att$basinSize)
      cat("  Basin size range:", min(basinSizes), "to", max(basinSizes), "\n")
      cat("  Total basin coverage:", sum(basinSizes), "\n")
    }
    
    # Check for reasonable distribution
    if (nAttractors == length(attractorSizes)) {
      cat("  Samples per attractor: SUSPICIOUS - every sample became unique attractor\n")
    } else {
      cat("  Sample convergence: GOOD - multiple samples converged to same attractors\n")
    }
  }
  
  cat("=====================================\n\n")
  
  return(list(
    network_genes = length(boolnet$genes),
    n_attractors = nAttractors,
    attractor_sizes = if(nAttractors > 0) attractorSizes else NULL,
    is_suspicious = nAttractors == length(attractorSizes)
  ))
}

#' Create attractor size distribution plot
#'
#' @param attractors Attractors object from getAttractors()
#' @return ggplot object
plotAttractorSizeDistribution <- function(attractors) {
  
  if (length(attractors$attractors) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No attractors found", size = 6) +
             theme_void())
  }
  
  attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))
  
  plotData <- data.frame(
    AttractorID = 1:length(attractorSizes),
    Size = attractorSizes,
    Type = ifelse(attractorSizes == 1, "Steady State", "Cycle")
  )
  
  p <- ggplot(plotData, aes(x = Size, fill = Type)) +
    geom_histogram(binwidth = 1, alpha = 0.7, color = "white") +
    scale_fill_manual(values = c("Steady State" = "steelblue", "Cycle" = "orange")) +
    theme_minimal() +
    labs(
      title = "Attractor Size Distribution",
      subtitle = paste("Total:", length(attractorSizes), "attractors"),
      x = "Attractor Size (Number of States)",
      y = "Count",
      fill = "Attractor Type"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

#' Create basin size distribution plot
#'
#' @param attractors Attractors object from getAttractors()
#' @return ggplot object
plotBasinSizeDistribution <- function(attractors) {
  
  if (length(attractors$attractors) == 0 || 
      !"basinSize" %in% names(attractors$attractors[[1]])) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "Basin size data not available", size = 6) +
             theme_void())
  }
  
  basinSizes <- sapply(attractors$attractors, function(att) att$basinSize)
  attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))
  
  plotData <- data.frame(
    AttractorID = 1:length(basinSizes),
    BasinSize = basinSizes,
    AttractorSize = attractorSizes,
    Type = ifelse(attractorSizes == 1, "Steady State", "Cycle")
  )
  
  # Sort by basin size for better visualization
  plotData <- plotData[order(-plotData$BasinSize), ]
  plotData$AttractorID <- factor(plotData$AttractorID, levels = plotData$AttractorID)
  
  p <- ggplot(plotData, aes(x = AttractorID, y = BasinSize, fill = Type)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("Steady State" = "steelblue", "Cycle" = "orange")) +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(
      title = "Basin of Attraction Sizes",
      subtitle = paste("Total basin coverage:", sum(basinSizes)),
      x = "Attractors (ordered by basin size)",
      y = "Basin Size",
      fill = "Attractor Type"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

#' Generate comprehensive attractor summary table
#'
#' @param attractors Attractors object from getAttractors()
#' @param boolnet BoolNet network object
#' @return Data frame with attractor details
generateAttractorSummaryTable <- function(attractors, boolnet) {
  
  if (length(attractors$attractors) == 0) {
    return(data.frame(Message = "No attractors found"))
  }
  
  summaryTable <- data.frame(
    AttractorID = integer(),
    AttractorSize = integer(),
    BasinSize = integer(),
    Type = character(),
    Stability = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(attractors$attractors)) {
    att <- attractors$attractors[[i]]
    
    attSize <- ncol(att$involvedStates)
    basinSize <- if("basinSize" %in% names(att)) att$basinSize else NA
    
    summaryTable <- rbind(summaryTable, data.frame(
      AttractorID = i,
      AttractorSize = attSize,
      BasinSize = basinSize,
      Type = ifelse(attSize == 1, "Steady State", "Cycle"),
      Stability = case_when(
        is.na(basinSize) ~ "Unknown",
        basinSize > 100 ~ "High",
        basinSize > 10 ~ "Medium", 
        TRUE ~ "Low"
      )
    ))
  }
  
  # Sort by basin size (largest first)
  summaryTable <- summaryTable[order(-summaryTable$BasinSize, na.last = TRUE), ]
  
  return(summaryTable)
}

#' Full validation report with plots and tables
#'
#' @param attractors Attractors object from getAttractors()
#' @param boolnet BoolNet network object
#' @param savePath Path to save plots (optional)
#' @return List with validation results
generateValidationReport <- function(attractors, boolnet, savePath = NULL) {
  
  # Run validation checks
  validation <- validateBoolNetResults(attractors, boolnet)
  
  # Generate summary table
  summaryTable <- generateAttractorSummaryTable(attractors, boolnet)
  
  # Create plots
  p1 <- plotAttractorSizeDistribution(attractors)
  p2 <- plotBasinSizeDistribution(attractors)
  
  # Save plots if path provided
  if (!is.null(savePath)) {
    ggsave(paste0(savePath, "attractor_sizes.png"), p1, width = 10, height = 6, dpi = 300)
    ggsave(paste0(savePath, "basin_sizes.png"), p2, width = 10, height = 6, dpi = 300)
    
    # Save summary table
    write.table(summaryTable, paste0(savePath, "attractor_summary.tsv"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("Validation plots and table saved to:", savePath, "\n")
  }
  
  # Print key findings
  cat("=== KEY FINDINGS ===\n")
  if (validation$is_suspicious) {
    cat("⚠️  SUSPICIOUS: Every sample became a unique attractor\n")
    cat("   This suggests the network has trivial dynamics\n")
  } else {
    cat("✅ GOOD: Samples converged to shared attractors\n")
    cat("   This suggests meaningful network dynamics\n")
  }
  
  if (validation$n_attractors > 0) {
    steadyStates <- sum(validation$attractor_sizes == 1)
    cycles <- sum(validation$attractor_sizes > 1)
    cat("📊 SUMMARY:", steadyStates, "steady states,", cycles, "cycles\n")
  }
  cat("====================\n")
  
  return(list(
    validation = validation,
    summary_table = summaryTable,
    plots = list(attractor_sizes = p1, basin_sizes = p2)
  ))
}

# =============================================================================
# Usage in Script 09 (add after attractor computation)
# =============================================================================

# # Add this validation section after successful attractor computation:
# 
# # --- Validate Results ---
# message("Validating BoolNet analysis results...")
# 
# if (config$saveResults) {
#   validationPath <- paste0(paths$base$plots, cellType, "_", trajectory, "_validation_")
#   
#   validationReport <- generateValidationReport(
#     attractors = attractors,
#     boolnet = boolnet, 
#     savePath = validationPath
#   )
#   
#   # Print summary table
#   cat("ATTRACTOR SUMMARY TABLE:\n")
#   print(head(validationReport$summary_table, 10))
#   
# } else {
#   # Just run basic validation without saving
#   validation <- validateBoolNetResults(attractors, boolnet)
# }