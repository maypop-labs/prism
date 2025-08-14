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
if (!file.exists(paths$static$mergedSeurat)) stop("Merged Seurat object not found")
if (config$verbose) { message("Loading merged Seurat RDS file") }
mergedSeurat <- readRDS(paths$static$mergedSeurat)

# --- Assign Cell Types Using SingleR ---
if (config$verbose) { message("Assigning cell types using SingleR") }
sce       <- as.SingleCellExperiment(mergedSeurat)
ref       <- celldex::HumanPrimaryCellAtlasData()
cellAnno  <- SingleR(test = sce, ref = ref, labels = ref$label.main)
mergedSeurat$cellType <- cellAnno$labels

# --- Sort By Most Common Cell Types ---
cellTypeFreq <- sort(table(mergedSeurat$cellType), decreasing = TRUE)
cellTypes <- names(cellTypeFreq[cellTypeFreq >= config$singleRMinimumNumberOfCells])
if (config$verbose) { message("Cell types: ", paste(cellTypes, collapse = ", ")) }

if (config$saveResults) {
  saveRDS(cellTypes, file = paste0(paths$base$rds, "cell_types.rds"))
}

# --- Process Each Cell Type ---
for (ct in cellTypes) {
  readableCt <- toTitleCase(gsub("_", " ", ct))
  message("Processing cell type: ", readableCt)

  ctObj <- subset(mergedSeurat, subset = cellType == ct)
  ctObj <- FindVariableFeatures(ctObj, nfeatures = config$singleRNumberOfFeatures)
  varGenes <- VariableFeatures(ctObj)
  ctObj <- ScaleData(ctObj, features = varGenes)
  ctObj <- RunPCA(ctObj, features = varGenes)
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

    if (config$verbose) { message("Saving top genes from PC1 and PC2 for ", readableCt) }
    loadings <- Loadings(ctObj, reduction = "pca")
    for (pc in 1:2) {
      pcName <- paste0("PC_", pc)
      topGenes <- names(sort(loadings[, pcName], decreasing = TRUE))[1:50]
      write.table(
        topGenes,
        file = paste0(paths$base$txt, ct, "_PCA_top_genes_", pcName, ".txt"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
    
    if (config$verbose) { message("Saving plots for ", readableCt) }
    
    ggsave(paste0(paths$base$plots, "figure3_", readableCt, ".png"), pPCA,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    ggsave(paste0(paths$base$plots, "figure4_", readableCt, ".png"), pUMAP,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    
    if (config$verbose) { message("Saving Seurat file for ", readableCt) }
    saveRDS(ctObj, file = paste0(paths$base$rds, ct, "_seurat.rds"))
  }
}

message("Done!")
