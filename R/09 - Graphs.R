# =============================================================================
# 09 - Graphs.R
# Purpose: Generate publication-ready visualizations of PRISM analysis results
# =============================================================================

# Load required packages and initialize paths
source("managers/attractorManager.R")
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")
config <- initializeScript()
pathInfo <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths <- pathInfo$paths
cellType <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths <- getCellTypeFilePaths(paths$base, cellType)
ptPaths <- getTrajectoryFilePaths(paths$base, cellType, trajectory)

# Ensure plots directory exists
if (!dir.exists(paths$base$plots)) {
  dir.create(paths$base$plots, recursive = TRUE)
}

message("Loading data for visualization...")

# Load perturbation analysis results (optional - script continues if missing)
singleKD <- loadObject(ptPaths$singleTargetsKD, config, "single gene KD results", required = FALSE)
singleOE <- loadObject(ptPaths$singleTargetsOE, config, "single gene OE results", required = FALSE)

# Load double perturbation results (optional)
doubleKD <- loadObject(ptPaths$doubleTargetsKD, config, "double KD results", required = FALSE)
doubleOE <- loadObject(ptPaths$doubleTargetsOE, config, "double OE results", required = FALSE)
doubleMix <- loadObject(ptPaths$doubleTargetsMix, config, "mixed KD/OE results", required = FALSE)

# Load reference data (optional)
networkSummary <- loadObject(ptPaths$networkAgingSummary, config, "network aging summary", required = FALSE)
initialScore <- if (!is.null(networkSummary)) networkSummary$OverallAgingScore else NULL

message("Creating visualizations from available data...")

# Initialize plot collections
attractorPlots <- list()
singleGenePlots <- list()
doubleGenePlots <- list()

# =============================================================================
# ATTRACTOR LANDSCAPE VISUALIZATION
# =============================================================================

# Load attractor analysis data (optional)
attractorDfScores <- loadObject(ptPaths$attractorScoresCombined, config, "combined attractor scores", required = FALSE)

if (!is.null(attractorDfScores) && nrow(attractorDfScores) > 0) {
  message("Creating attractor landscape visualization...")
  
  # Add attractor index for labeling
  attractorDfScores$AttractorIndex <- seq_len(nrow(attractorDfScores))
  
  # Calculate plot limits with padding
  ageRange <- range(attractorDfScores$AgeScore, na.rm = TRUE)
  agePadding <- diff(ageRange) * 0.1  # 10% padding
  if (agePadding < 0.01) agePadding <- 0.05  # Minimum padding for single points
  
  basinRange <- range(attractorDfScores$BasinSize, na.rm = TRUE)
  basinPadding <- diff(basinRange) * 0.1  # 10% padding
  if (basinPadding < 0.01) basinPadding <- 0.05  # Minimum padding for single points
  
  # Determine if we need special handling for single attractor
  numAttractors <- nrow(attractorDfScores)
  
  # Enhanced bubble plot showing three dimensions of attractor properties
  if (numAttractors == 1) {
    # Special handling for single attractor - use fixed large size
    attractorLandscapePlot <- ggplot(attractorDfScores, aes(x = AgeScore, y = BasinSize, fill = AttractorScore)) +
      # Single large point with black border
      geom_point(size = 30, alpha = 1.0, shape = 21, color = "black", stroke = 3) +
      # Text label showing attractor index
      geom_text(aes(label = AttractorIndex), 
                color = "white", 
                size = 8, 
                fontface = "bold",
                show.legend = FALSE) +

      scale_x_continuous(
        name = "Age Score (0 = Young-like, 1 = Old-like)",
        limits = c(ageRange[1] - agePadding, ageRange[2] + agePadding),
        expand = c(0.02, 0)
      ) +
      scale_y_continuous(
        name = "Basin Size Fraction",
        limits = c(max(0, basinRange[1] - basinPadding), basinRange[2] + basinPadding),
        expand = c(0.02, 0)
      ) +
      scale_fill_gradient2(
        name = "Composite\nAging Score",
        low = "#0571b0",      # Strong blue (young-like)
        mid = "#f4a582",      # Visible peach (neutral)
        high = "#ca0020",     # Strong red (old-like)
        midpoint = median(attractorDfScores$AttractorScore, na.rm = TRUE),
        limits = c(0, max(5, attractorDfScores$AttractorScore * 1.2)),
        guide = guide_colorbar(title.position = "top")
      ) +
      labs(
        title = "Attractor Landscape: Aging Bias vs. Network Dominance",
        subtitle = "Single attractor dominates entire state space",
        x = NULL,
        y = NULL,
        caption = "Point size fixed for visibility | Color = Composite aging contribution | Number = Attractor index"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 9, color = "gray50"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  } else {
    # Multiple attractors - use size aesthetic
    attractorLandscapePlot <- ggplot(attractorDfScores, aes(x = AgeScore, y = BasinSize, 
                                                            size = Stability, fill = AttractorScore)) +
      geom_point(alpha = 1.0, shape = 21, color = "black", stroke = 2) +
      geom_text(aes(label = AttractorIndex), 
                color = "white", 
                size = 5, 
                fontface = "bold",
                show.legend = FALSE) +
      scale_x_continuous(
        name = "Age Score (0 = Young-like, 1 = Old-like)",
        limits = c(ageRange[1] - agePadding, ageRange[2] + agePadding),
        expand = c(0.02, 0)
      ) +
      scale_y_continuous(
        name = "Basin Size Fraction",
        limits = c(max(0, basinRange[1] - basinPadding), basinRange[2] + basinPadding),
        expand = c(0.02, 0)
      ) +
      scale_size_continuous(
        name = "Stability",
        range = c(10, 25),
        guide = "none"
      ) +
      scale_fill_gradient2(
        name = "Composite\nAging Score",
        low = "#0571b0",
        mid = "#f4a582",
        high = "#ca0020",
        midpoint = median(attractorDfScores$AttractorScore, na.rm = TRUE),
        guide = guide_colorbar(title.position = "top")
      ) +
      labs(
        title = "Attractor Landscape: Aging Bias vs. Network Dominance",
        x = NULL,
        y = NULL,
        caption = "Point size = Stability | Color = Composite aging contribution | Numbers = Attractor index"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 9, color = "gray50"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  }
  
  attractorPlots$attractor_landscape <- attractorLandscapePlot
  message("Attractor landscape visualization created (", nrow(attractorDfScores), " attractors)")
} else {
  message("Skipping attractor landscape visualization (data not available)")
}

# =============================================================================
# GENESWITCHES TIMELINE VISUALIZATION
# =============================================================================

# Load GeneSwitches results (optional)
geneSwitchesData <- loadObject(ptPaths$geneSwitchesTsv, config, "GeneSwitches results", required = FALSE)

if (!is.null(geneSwitchesData) && nrow(geneSwitchesData) > 0) {
  message("Creating GeneSwitches timeline visualization...")
  
  # Filter genes by quality metrics
  # - FDR < 0.05
  # - Pseudo R² > 0.03
  filteredSwitches <- geneSwitchesData[
    geneSwitchesData$fdr < 0.05 & 
    geneSwitchesData$pseudoR2s > 0.03,
  ]
  
  if (nrow(filteredSwitches) == 0) {
    message("No genes passed GeneSwitches filtering criteria (FDR < 0.05, Pseudo R² > 0.03). Skipping timeline plot.")
  } else {
    # Select top 30 genes by Pseudo R² (quality of fit)
    filteredSwitches <- filteredSwitches[order(-filteredSwitches$pseudoR2s), ]
    top30Switches <- head(filteredSwitches, 30)
    
    # Separate upregulated and downregulated genes for positioning
    top30Switches$PlotPosition <- ifelse(
      top30Switches$direction == "up",
      top30Switches$pseudoR2s,    # Above timeline (positive y)
      -top30Switches$pseudoR2s    # Below timeline (negative y)
    )
    
    # Create the timeline plot
    geneSwitchesTimelinePlot <- ggplot(top30Switches, 
                                        aes(x = pseudotime, y = PlotPosition, label = geneId)) +
      # Horizontal timeline at y = 0
      geom_hline(yintercept = 0, color = "black", linewidth = 1) +
      
      # Gene points
      geom_point(aes(color = direction), size = 3, alpha = 0.8) +
      
      # Gene labels
      ggrepel::geom_text_repel(aes(color = direction), size = 3, 
                               max.overlaps = Inf, 
                               segment.size = 0.2,
                               box.padding = 0.5,
                               point.padding = 0.3,
                               show.legend = FALSE) +
      
      # Color scheme: up = green, down = red
      scale_color_manual(
        name = "Switch Direction",
        values = c("up" = "#1a9850", "down" = "#d73027"),
        labels = c("up" = "Upregulated", "down" = "Downregulated")
      ) +
      
      # Axis labels and title
      labs(
        title = paste0("Gene Expression Switches Along Pseudotime"),
        subtitle = paste0(cellType, " - ", trajectory, " (Top 30 by Pseudo R²)"),
        x = "Pseudotime",
        y = "Quality of Fitting (Pseudo R²)",
        caption = "Genes above timeline are upregulated; genes below are downregulated"
      ) +
      
      # Theme
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5)
      )
    
    attractorPlots$geneSwitches <- geneSwitchesTimelinePlot
    message("GeneSwitches timeline visualization created (", nrow(top30Switches), " genes)")
  }
} else {
  message("Skipping GeneSwitches timeline visualization (data not available)")
}

# =============================================================================
# PERTURBATION RESULTS PLOTS
# =============================================================================

# Publication theme
applyPublicationTheme <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}

# Single Gene Perturbation Results
if (!is.null(singleKD) && !is.null(singleOE) && nrow(singleKD) > 0 && nrow(singleOE) > 0) {
  message("Creating single gene perturbation visualizations...")
  
  # Combine single gene results
  singleKD$Mode <- "Knockdown"
  singleOE$Mode <- "Overexpression"
  singleCombined <- rbind(singleKD, singleOE)
  
  # Filter for successful results and add improvement categories
  successfulSingle <- singleCombined[singleCombined$Success & !is.na(singleCombined$Delta), ]
  
  if (nrow(successfulSingle) > 0) {
    successfulSingle$Improvement <- ifelse(successfulSingle$Delta < -0.01, "Significant", 
                                    ifelse(successfulSingle$Delta < 0, "Moderate", "None"))
    successfulSingle$Improvement <- factor(successfulSingle$Improvement, 
                                    levels = c("Significant", "Moderate", "None"))
    
    # Plot: Top 20 single gene effects
    numToShow <- min(20, nrow(successfulSingle))
    top20Single <- head(successfulSingle[order(successfulSingle$Delta), ], numToShow)
    
    singleEffectsPlot <- ggplot(top20Single, aes(x = reorder(Gene, Delta), y = Delta, fill = Mode)) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_hline(yintercept = -0.01, linetype = "dotted", color = "red", alpha = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Knockdown" = "#d73027", "Overexpression" = "#1a9850")) +
      labs(
        title = paste0("Top ", numToShow, " Single Gene Perturbation Effects"),
        x = "Gene",
        y = "Change in Aging Score (Delta)",
        fill = "Intervention",
        caption = "Negative values indicate aging reduction (improvement)"
      ) +
      applyPublicationTheme()
    
    singleGenePlots$single_effects <- singleEffectsPlot
    message("Single gene perturbation visualization created (", numToShow, " genes)")
  } else {
    message("No successful single gene perturbations found. Skipping single gene plots.")
  }
} else {
  message("Skipping single gene perturbation visualizations (data not available)")
}

# Double Gene Perturbation Results (if available)
if (!is.null(doubleKD) || !is.null(doubleOE) || !is.null(doubleMix)) {
  message("Creating double gene perturbation visualizations...")
  
  # Combine double perturbation results
  doubleResults <- list()
  
  if (!is.null(doubleKD) && nrow(doubleKD) > 0) {
    doubleKD$Type <- "Double KD"
    doubleKD$GeneCombo <- paste(doubleKD$Gene1, "+", doubleKD$Gene2)
    # Standardize columns for mixing
    doubleKD$KD_Gene <- doubleKD$Gene1
    doubleKD$OE_Gene <- doubleKD$Gene2
    doubleResults <- append(doubleResults, list(doubleKD))
  }
  
  if (!is.null(doubleOE) && nrow(doubleOE) > 0) {
    doubleOE$Type <- "Double OE"
    doubleOE$GeneCombo <- paste(doubleOE$Gene1, "+", doubleOE$Gene2)
    # Standardize columns for mixing
    doubleOE$KD_Gene <- doubleOE$Gene1
    doubleOE$OE_Gene <- doubleOE$Gene2
    doubleResults <- append(doubleResults, list(doubleOE))
  }
  
  if (!is.null(doubleMix) && nrow(doubleMix) > 0) {
    doubleMix$Type <- "Mixed KD/OE"
    doubleMix$GeneCombo <- paste(doubleMix$KD_Gene, "(KD) +", doubleMix$OE_Gene, "(OE)")
    doubleMix$Gene1 <- doubleMix$KD_Gene
    doubleMix$Gene2 <- doubleMix$OE_Gene
    doubleResults <- append(doubleResults, list(doubleMix))
  }
  
  if (length(doubleResults) > 0) {
    # Ensure all data frames have the same columns before combining
    commonCols <- c("AgingScore", "Delta", "Type", "GeneCombo")
    
    # Add missing columns to each data frame
    for (i in seq_along(doubleResults)) {
      for (col in commonCols) {
        if (!col %in% colnames(doubleResults[[i]])) {
          doubleResults[[i]][[col]] <- NA
        }
      }
      # Keep only the common columns
      doubleResults[[i]] <- doubleResults[[i]][, commonCols, drop = FALSE]
    }
    
    doubleCombined <- do.call(rbind, doubleResults)
    doubleCombined$Improvement <- ifelse(doubleCombined$Delta < -0.01, "Significant", 
                                  ifelse(doubleCombined$Delta < 0, "Moderate", "None"))
    
    # Plot: Top double gene combinations
    numToShow <- min(15, nrow(doubleCombined))
    top15Double <- head(doubleCombined[order(doubleCombined$Delta), ], numToShow)
    
    doubleEffectsPlot <- ggplot(top15Double, aes(x = reorder(GeneCombo, Delta), y = Delta, fill = Type)) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_hline(yintercept = -0.01, linetype = "dotted", color = "red", alpha = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Double KD" = "#d73027", "Double OE" = "#1a9850", "Mixed KD/OE" = "#4575b4")) +
      labs(
        title = paste0("Top ", numToShow, " Double Gene Perturbation Effects"),
        x = "Gene Combination",
        y = "Change in Aging Score (Delta)",
        fill = "Intervention Type",
        caption = "Negative values indicate aging reduction (improvement)"
      ) +
      applyPublicationTheme() +
      theme(axis.text.y = element_text(size = 8))
    
    doubleGenePlots$double_effects <- doubleEffectsPlot
    message("Double gene perturbation visualization created (", numToShow, " combinations)")
  }
} else {
  message("Skipping double gene perturbation visualizations (data not available)")
}

# Combine all plots
allVisualizationPlots <- c(attractorPlots, singleGenePlots, doubleGenePlots)

if (length(allVisualizationPlots) > 0) {
  message("All visualizations created: ", length(allVisualizationPlots), " plots")
} else {
  message("No visualizations could be created - all required data files are missing")
}

# =============================================================================
# SAVE PLOTS
# =============================================================================

if (config$saveResults && length(allVisualizationPlots) > 0) {
  message("Saving visualization plots...")
  
  savedPlots <- character()
  for (plotName in names(allVisualizationPlots)) {
    plotObj <- allVisualizationPlots[[plotName]]
    
    if (!is.null(plotObj)) {
      fileName <- paste0(cellType, "_", trajectory, "_", plotName, ".png")
      filePath <- file.path(paths$base$plots, fileName)
      
      # Use appropriate dimensions based on plot type
      plotWidth <- if (grepl("effects", plotName)) 10 else 8
      plotHeight <- if (grepl("effects", plotName)) 8 else 6
      
      ggsave(filePath, plot = plotObj, width = plotWidth, height = plotHeight, 
             dpi = config$figDPI %||% 300, bg = "white")
      savedPlots <- c(savedPlots, fileName)
    }
  }
  
  if (length(savedPlots) > 0) {
    message("Saved ", length(savedPlots), " visualization plots:")
    for (plot in savedPlots) {
      message("  - ", plot)
    }
  } else {
    message("No plots were saved")
  }
} else {
  if (!config$saveResults) {
    message("Skipping plot saving (saveResults = FALSE)")
  } else {
    message("No plots to save")
  }
}

message("Done!")
