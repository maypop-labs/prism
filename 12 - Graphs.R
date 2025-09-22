# =============================================================================
# 12 - Graphs.R
# Purpose: Generate visualizations of PRISM analysis results
# =============================================================================

# Load required packages and initialize paths
source("managers/attractorManager.R")
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")
config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)

# Ensure plots directory exists
if (!dir.exists(paths$base$plots)) {
  dir.create(paths$base$plots, recursive = TRUE)
}

message("Loading data for visualization...")

# Load data objects
seuratMerged <- loadObject(paths$static$seuratMerged, config, "merged Seurat object")
cds <- loadMonocle3(ptPaths$monocle3, config, "pseudotime trajectory")

# Create generic gene mapping for IP protection
allGenes <- rownames(seuratMerged)
set.seed(config$randomSeed)  # Reproducible mapping
genericNames <- paste0("Gene", sprintf("%04d", sample(1:9999, length(allGenes))))
geneMapping <- setNames(genericNames, allGenes)

message("Creating visualizations...")

# Publication theme
applyPublicationTheme <- function(plot) {
  plot + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# 1. UMAP Plot by Age
umapPlot <- NULL
if ("umap" %in% names(seuratMerged@reductions)) {
  umapPlot <- DimPlot(seuratMerged, reduction = "umap", group.by = "age", pt.size = 0.5) +
    ggtitle("UMAP Embedding by Donor Age") +
    labs(x = "UMAP1", y = "UMAP2", color = "Age") +
    applyPublicationTheme()
} else if ("umap.unintegrated" %in% names(seuratMerged@reductions)) {
  umapPlot <- DimPlot(seuratMerged, reduction = "umap.unintegrated", group.by = "age", pt.size = 0.5) +
    ggtitle("UMAP Embedding by Donor Age") +
    labs(x = "UMAP1", y = "UMAP2", color = "Age") +
    applyPublicationTheme()
}

# 2. Pseudotime Trajectory
pseudotimePlot <- plot_cells(cds,
         color_cells_by = "pseudotime",
         label_cell_groups = FALSE,
         show_trajectory_graph = TRUE) +
  ggtitle("Pseudotime Trajectory Analysis") +
  applyPublicationTheme()

# 3. Cell Count by Age
cellCounts <- table(seuratMerged$age)
cellCountDf <- data.frame(
  Age = names(cellCounts),
  CellCount = as.numeric(cellCounts)
)

cellCountPlot <- ggplot(cellCountDf, aes(x = Age, y = CellCount, fill = Age)) +
  geom_col(alpha = 0.8) +
  scale_fill_viridis_d() +
  ggtitle("Cell Count Distribution by Age") +
  labs(x = "Donor Age", y = "Number of Cells") +
  applyPublicationTheme() +
  theme(legend.position = "none")

# Save plots
if (config$saveResults) {
  plotList <- list(
    umap_plot = umapPlot,
    pseudotime_plot = pseudotimePlot,
    cell_count_plot = cellCountPlot
  )
  
  # Remove NULL plots
  plotList <- plotList[!sapply(plotList, is.null)]
  
  savedPlots <- character()
  for (plotName in names(plotList)) {
    plotObj <- plotList[[plotName]]
    
    if (!is.null(plotObj)) {
      fileName <- paste0(cellType, "_", trajectory, "_", plotName, ".png")
      filePath <- file.path(paths$base$plots, fileName)
      
      ggsave(filePath, plot = plotObj, width = 10, height = 8, dpi = 300)
      savedPlots <- c(savedPlots, fileName)
    }
  }
  
}

message("Done!")
