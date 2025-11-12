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

# Load perturbation analysis results
singleKD <- loadObject(ptPaths$singleTargetsKD, config, "single gene KD results")
singleOE <- loadObject(ptPaths$singleTargetsOE, config, "single gene OE results")

# Load double perturbation results (if they exist)
doubleKD <- loadObject(ptPaths$doubleTargetsKD, config, "double KD results", required = FALSE)
doubleOE <- loadObject(ptPaths$doubleTargetsOE, config, "double OE results", required = FALSE)
doubleMix <- loadObject(ptPaths$doubleTargetsMix, config, "mixed KD/OE results", required = FALSE)

# Load reference data for comparison
networkSummary <- loadObject(ptPaths$networkAgingSummary, config, "network aging summary")
initialScore <- networkSummary$OverallAgingScore

message("Creating perturbation results visualizations...")

# =============================================================================
# ATTRACTOR LANDSCAPE VISUALIZATION
# =============================================================================

# Load attractor analysis data
attractorDfScores <- loadObject(ptPaths$attractorScoresCombined, config, "combined attractor scores")

# Enhanced bubble plot showing three dimensions of attractor properties
attractorLandscapePlot <- ggplot(attractorDfScores, aes(x = AgeScore, y = BasinSize, 
                                                        size = Stability, color = AttractorScore)) +
  geom_point(alpha = 1.0, stroke = 1, color = "black") +
  geom_point(alpha = 1.0, stroke = 0) +
  scale_y_continuous(
    name = "Basin Size Fraction"
  ) +
  scale_size_continuous(
    name = "Stability",
    range = c(5, 15),
    guide = "none"
  ) +
  scale_color_gradient2(
    name = "Composite\nAging Score",
    low = "#2166ac",
    mid = "#f7f7f7", 
    high = "#b2182b",
    midpoint = median(attractorDfScores$AttractorScore, na.rm = TRUE),
    guide = guide_colorbar(title.position = "top")
  ) +
  labs(
    title = "Attractor Landscape: Aging Bias vs. Network Dominance",
    x = "Age Score (0 = Young-like, 1 = Old-like)",
    y = "Basin Size Fraction",
    caption = "Point size = Stability | Color = Composite aging contribution"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 9, color = "gray50")
  )

# Create attractor plot list
attractorPlots <- list(
  attractor_landscape = attractorLandscapePlot
)

message("Attractor landscape visualization created")

# =============================================================================
# GENESWITCHES TIMELINE VISUALIZATION
# =============================================================================

message("Creating GeneSwitches timeline visualization...")

# Load GeneSwitches results
geneSwitchesData <- loadObject(ptPaths$geneSwitchesTsv, config, "GeneSwitches results")

# Filter genes by quality metrics
# - FDR < 0.05
# - Pseudo R² > 0.03
filteredSwitches <- geneSwitchesData[
  geneSwitchesData$fdr < 0.05 & 
  geneSwitchesData$pseudoR2s > 0.03,
]

if (nrow(filteredSwitches) == 0) {
  warning("No genes passed GeneSwitches filtering criteria. Skipping timeline plot.")
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
  
  # Add to plot list
  attractorPlots$geneSwitches <- geneSwitchesTimelinePlot
  
  message("GeneSwitches timeline visualization created (", nrow(top30Switches), " genes)")
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

# 1. Single Gene Perturbation Results
# Combine single gene results
singleKD$Mode <- "Knockdown"
singleOE$Mode <- "Overexpression"
singleCombined <- rbind(singleKD, singleOE)

# Filter for successful results and add improvement categories
successfulSingle <- singleCombined[singleCombined$Success & !is.na(singleCombined$Delta), ]
successfulSingle$Improvement <- ifelse(successfulSingle$Delta < -0.01, "Significant", 
                                ifelse(successfulSingle$Delta < 0, "Moderate", "None"))
successfulSingle$Improvement <- factor(successfulSingle$Improvement, 
                                levels = c("Significant", "Moderate", "None"))

# Plot 1: Top 20 single gene effects
top20Single <- head(successfulSingle[order(successfulSingle$Delta), ], 20)

singleEffectsPlot <- ggplot(top20Single, aes(x = reorder(Gene, Delta), y = Delta, fill = Mode)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -0.01, linetype = "dotted", color = "red", alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Knockdown" = "#d73027", "Overexpression" = "#1a9850")) +
  labs(
    title = "Top 20 Single Gene Perturbation Effects",
    x = "Gene",
    y = "Change in Aging Score (Delta)",
    fill = "Intervention",
    caption = "Negative values indicate aging reduction (improvement)"
  ) +
  applyPublicationTheme()

# Create plot list for single gene results
singleGenePlots <- list(
  single_effects = singleEffectsPlot
)

# 4. Double Gene Perturbation Results (if available)
doubleGenePlots <- list()

if (!is.null(doubleKD) && nrow(doubleKD) > 0) {
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
    
    # Plot 4: Top double gene combinations
    top15Double <- head(doubleCombined[order(doubleCombined$Delta), ], 15)
    
    doubleEffectsPlot <- ggplot(top15Double, aes(x = reorder(GeneCombo, Delta), y = Delta, fill = Type)) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_hline(yintercept = -0.01, linetype = "dotted", color = "red", alpha = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Double KD" = "#d73027", "Double OE" = "#1a9850", "Mixed KD/OE" = "#4575b4")) +
      labs(
        title = "Top 15 Double Gene Perturbation Effects",
        x = "Gene Combination",
        y = "Change in Aging Score (Delta)",
        fill = "Intervention Type",
        caption = "Negative values indicate aging reduction (improvement)"
      ) +
      applyPublicationTheme() +
      theme(axis.text.y = element_text(size = 8))
    
    doubleGenePlots <- list(
      double_effects = doubleEffectsPlot
    )
  }
}

# Combine all plots
allVisualizationPlots <- c(attractorPlots, singleGenePlots, doubleGenePlots)

message("All visualizations created: ", length(allVisualizationPlots), " plots")

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
  
  message("Saved ", length(savedPlots), " visualization plots:")
  for (plot in savedPlots) {
    message("  - ", plot)
  }
}

message("Done!")
