# =============================================================================
# 02 - Cell Type Identification
#
# Load merged Seurat object, assign cell types using SingleR and
# HumanPrimaryCellAtlas, normalize donor representation within each cell type,
# and generate PCA/UMAP plots for each cell type.
# =============================================================================

# --- Initialization ---
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")
source("managers/normalizationManager.R")
config <- initializeScript()
pathInfo <- initializeInteractivePaths()
paths <- pathInfo$paths
ensureProjectDirectories(paths)
clearConsole()

# --- Load Data ---
seuratMerged <- loadObject(paths$static$seuratMerged, config, "merged Seurat object")

# --- Assign Cell Types Using SingleR ---
if (config$verbose) { message("Assigning cell types using SingleR") }
sce <- as.SingleCellExperiment(seuratMerged)
ref <- celldex::HumanPrimaryCellAtlasData()
cellAnno <- SingleR(test = sce, ref = ref, labels = ref$label.main)
seuratMerged$cellType <- cellAnno$labels

# --- Sort By Most Common Cell Types ---
cellTypeFreq <- sort(table(seuratMerged$cellType), decreasing = TRUE)
cellTypes <- names(cellTypeFreq[cellTypeFreq >= config$singleRMinimumNumberOfCells])
if (config$verbose) { message("Cell types: ", paste(cellTypes, collapse = ", ")) }

# --- Normalize Donor Representation Within Cell Types ---
if (config$verbose) { 
  message("\n=== Pre-Normalization Donor Statistics ===") 
  print(reportDonorStats(seuratMerged, donorColumn = "donorID", cellTypeColumn = "cellType"))
}

seuratMerged <- normalizeCellTypeByDonor(
  seuratObj = seuratMerged,
  donorColumn = "donorID",
  cellTypeColumn = "cellType",
  floorCells = config$normFloorCells,
  capCells = config$normCapCells,
  quantileCut = config$normQuantileCut,
  numReplicates = 1,
  seedValue = config$randomSeed
)

if (config$verbose) { 
  message("\n=== Post-Normalization Donor Statistics ===") 
  print(reportDonorStats(seuratMerged, donorColumn = "donorID", cellTypeColumn = "cellType"))
}

# --- Process Each Cell Type ---
for (ct in cellTypes) {
  readableCt <- toTitleCase(gsub("_", " ", ct))
  if (config$verbose) { message("Processing cell type: ", readableCt) }
  ctPaths <- getCellTypeFilePaths(paths$base, ct)
  ctObj <- subset(seuratMerged, subset = cellType == ct)
  ctObj <- DietSeurat(ctObj, scale.data = FALSE)
  ctObj <- FindVariableFeatures(ctObj, nfeatures = config$singleRNumberOfFeatures)
  varFt <- VariableFeatures(ctObj)
  ctObj <- ScaleData(ctObj, features = varFt)
  ctObj <- RunPCA(ctObj, features = varFt)
  ctObj <- FindNeighbors(ctObj, dims = 1:10)
  ctObj <- FindClusters(ctObj, resolution = 1)
  ctObj <- RunUMAP(ctObj, dims = 1:10)

  # --- Plotting ---
  pPCA <- DimPlot(ctObj,
                  reduction = "pca",
                  group.by = "age",
                  pt.size = config$pointSize,
                  stroke.size = config$strokeSize,
                  alpha = config$plotAlpha) +
                  ggtitle(paste("PCA for", readableCt)) +
                  theme(plot.title = element_text(hjust = config$hjust)) +
                  labs(x = "PC1") +
                  labs(y = "PC2")

  pUMAP <- DimPlot(ctObj,
                   reduction = "umap",
                   group.by = "age",
                   pt.size = config$pointSize,
                   stroke.size = config$strokeSize,
                   alpha = config$plotAlpha) +
                   ggtitle(paste("UMAP for", readableCt)) +
                   theme(plot.title = element_text(hjust = config$hjust)) +
                   labs(x = "UMAP1") +
                   labs(y = "UMAP2")

  # --- Save Results ---
  if (config$saveResults) {
    
    if (config$verbose) { message("Saving plots for ", readableCt) }
    
    ggsave(ctPaths$pcaPlot, pPCA,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    ggsave(ctPaths$umapPlot, pUMAP,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    
    saveObject(ctObj, ctPaths$seuratObject, config, "cell type Seurat object")
  }
}

# --- Plotting ---
pPCA <- DimPlot(seuratMerged,
                reduction = "pca",
                group.by = "cellType",
                pt.size = config$pointSize,
                stroke.size = config$strokeSize,
                alpha = config$plotAlpha) +
                ggtitle("PCA Using All Cells") +
                theme(plot.title = element_text(hjust = config$hjust)) +
                labs(x = "PC1") +
                labs(y = "PC2")

pUMAP <- DimPlot(seuratMerged,
                 reduction = "umap.unintegrated",
                 group.by = "cellType",
                 pt.size = config$pointSize,
                 stroke.size = config$strokeSize,
                 alpha = config$plotAlpha) +
                 ggtitle("UMAP Using All Cells") +
                 theme(plot.title = element_text(hjust = config$hjust)) +
                 labs(x = "UMAP1") +
                 labs(y = "UMAP2")

# --- Save Results ---
if (config$saveResults) {
  saveObject(cellTypes, paths$static$cellTypes, config, "cell types")
  
  cellTypeFreqDf <- data.frame(
    cell_type = names(cellTypeFreq),
    number = as.integer(cellTypeFreq),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  saveObject(cellTypeFreqDf, paths$static$cellTypeFreq, config, "cell type frequency table")
  if (config$verbose) { message("Saving plots") }
  ggsave(paths$static$pcaAllCellsByType, pPCA,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
  ggsave(paths$static$umapAllCellsByType, pUMAP,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
}

message("Done!")
