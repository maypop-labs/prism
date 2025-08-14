# =============================================================================
# 02 - Cell Type Identification
#
# Load merged Seurat object, assign cell types using SingleR and
# HumanPrimaryCellAtlas, and generate PCA/UMAP plots for each cell type.
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

config   <- initializeScript()
pathInfo <- initializeInteractivePaths()
paths    <- pathInfo$paths
ensureProjectDirectories(paths)
clearConsole()

# --- Load Merged Seurat Object ---
if (!file.exists(paths$static$seuratMerged)) stop("Merged Seurat object not found")
if (config$verbose) { message("Loading merged Seurat RDS file") }
seuratMerged <- readRDS(paths$static$seuratMerged)

# --- Assign Cell Types Using SingleR ---
if (config$verbose) { message("Assigning cell types using SingleR") }
sce       <- as.SingleCellExperiment(seuratMerged)
ref       <- celldex::HumanPrimaryCellAtlasData()
cellAnno  <- SingleR(test = sce, ref = ref, labels = ref$label.main)
seuratMerged$cellType <- cellAnno$labels

# --- Sort By Most Common Cell Types ---
cellTypeFreq <- sort(table(seuratMerged$cellType), decreasing = TRUE)
cellTypes    <- names(cellTypeFreq[cellTypeFreq >= config$singleRMinimumNumberOfCells])
if (config$verbose) { message("Cell types: ", paste(cellTypes, collapse = ", ")) }

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
                  reduction   = "pca",
                  group.by    = "age",
                  pt.size     = config$pointSize,
                  stroke.size = config$strokeSize,
                  alpha       = config$plotAlpha) +
                  ggtitle(paste("PCA for", readableCt)) +
                  theme(plot.title = element_text(hjust = config$hjust)) +
                  labs(x = "PC1") +
                  labs(y = "PC2")

  pUMAP <- DimPlot(ctObj,
                   reduction   = "umap",
                   group.by    = "age",
                   pt.size     = config$pointSize,
                   stroke.size = config$strokeSize,
                   alpha       = config$plotAlpha) +
                   ggtitle(paste("UMAP for", readableCt)) +
                   theme(plot.title = element_text(hjust = config$hjust)) +
                   labs(x = "UMAP1") +
                   labs(y = "UMAP2")

  # --- Save Results ---
  if (config$saveResults) {
    
    if (config$verbose) { message("Saving plots for ", readableCt) }
    
    ggsave(ctPaths$pcaPlot, pPCA,
           width  = config$figWidth,
           height = config$figHeight,
           dpi    = config$figDPI,
           units  = "in")
    ggsave(ctPaths$umapPlot, pUMAP,
           width  = config$figWidth,
           height = config$figHeight,
           dpi    = config$figDPI,
           units  = "in")
    
    if (config$verbose) { message("Saving Seurat file for ", readableCt) }
    saveRDS(ctObj, file = ctPaths$seuratObject)
  }
}

# --- Plotting ---
pPCA <- DimPlot(seuratMerged,
                reduction   = "pca",
                group.by    = "cellType",
                pt.size     = config$pointSize,
                stroke.size = config$strokeSize,
                alpha       = config$plotAlpha) +
                ggtitle("PCA Using All Cells") +
                theme(plot.title = element_text(hjust = config$hjust)) +
                labs(x = "PC1") +
                labs(y = "PC2")

pUMAP <- DimPlot(seuratMerged,
                 reduction   = "umap.unintegrated",
                 group.by    = "cellType",
                 pt.size     = config$pointSize,
                 stroke.size = config$strokeSize,
                 alpha       = config$plotAlpha) +
                 ggtitle("UMAP Using All Cells") +
                 theme(plot.title = element_text(hjust = config$hjust)) +
                 labs(x = "UMAP1") +
                 labs(y = "UMAP2")

# --- Save Results ---
if (config$saveResults) {
  if (config$verbose) { message("Saving cell types RDS file") }
  saveRDS(cellTypes, file = paths$static$cellTypes)
  
  if (config$verbose) { message("Saving cell type frequency table as TSV file") }
  cellTypeFreqDf <- data.frame(
    cell_type = names(cellTypeFreq),
    number = as.integer(cellTypeFreq),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(
    cellTypeFreqDf,
    file = file.path(paths$static$cellTypeFreq),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (config$verbose) { message("Saving plots") }
  ggsave(paths$static$pcaAllCellsByType, pPCA,
         width  = config$figWidth,
         height = config$figHeight,
         dpi    = config$figDPI,
         units  = "in")
  ggsave(paths$static$umapAllCellsByType, pUMAP,
         width  = config$figWidth,
         height = config$figHeight,
         dpi    = config$figDPI,
         units  = "in")
}

message("Done!")
