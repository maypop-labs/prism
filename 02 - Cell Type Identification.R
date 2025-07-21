# =============================================================================
# 02 - Cell Type Identification
#
# Load merged Seurat object, assign cell types using SingleR and
# HumanPrimaryCellAtlas, and generate PCA/UMAP plots for each cell type.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(Seurat)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)
library(tools)

# --- Config and options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
plotPath       <- paste0(config$rootPath, "results/plots/")
rdsPath        <- paste0(config$rootPath, "results/rds/")
tsvPath        <- paste0(config$rootPath, "results/tsv/")
txtPath        <- paste0(config$rootPath, "results/txt/")
seuratFile     <- paste0(rdsPath, "merged_seurat.rds")

dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,  recursive = TRUE, showWarnings = FALSE)

# --- Load Merged Seurat Object ---
if (!file.exists(seuratFile)) stop("Merged Seurat object not found")
mergedSeurat <- readRDS(seuratFile)

cat("\014")
cat("\n")

# --- Assign Cell Types Using SingleR ---
message("Assigning cell types using SingleR")
sce       <- as.SingleCellExperiment(mergedSeurat)
ref       <- celldex::HumanPrimaryCellAtlasData()
cellAnno  <- SingleR(test = sce, ref = ref, labels = ref$label.main)
mergedSeurat$cellType <- cellAnno$labels

# --- Sort By Most Common Cell Types ---
cellTypeFreq <- sort(table(mergedSeurat$cellType), decreasing = TRUE)
cellTypes <- names(cellTypeFreq[cellTypeFreq >= config$minCellTypeNumber])
message("Cell types: ", paste(cellTypes, collapse = ", "))

if (config$saveResults) {
  saveRDS(cellTypes, file = paste0(rdsPath, "cell_types.rds"))
}

# --- Process Each Cell Type ---
for (ct in cellTypes) {
  readableCt <- toTitleCase(gsub("_", " ", ct))
  message("=== Processing cell type: ", readableCt, " ===")

  ctObj <- subset(mergedSeurat, subset = cellType == ct)
  ctObj <- FindVariableFeatures(ctObj, nfeatures = 2000)
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
    
    # save cellTypeFreq table as tsv file
    cellTypeFreqDf <- data.frame(
      cell_type = names(cellTypeFreq),
      number = as.integer(cellTypeFreq),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    write.table(
      cellTypeFreqDf,
      file = file.path(tsvPath, "cell_types.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    message("Saving top genes from PC1 and PC2 for ", readableCt)
    loadings <- Loadings(ctObj, reduction = "pca")
    for (pc in 1:2) {
      pcName <- paste0("PC_", pc)
      topGenes <- names(sort(loadings[, pcName], decreasing = TRUE))[1:50]
      write.table(
        topGenes,
        file = paste0(txtPath, ct, "_PCA_top_genes_", pcName, ".txt"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
    
    message("Saving plots for ", readableCt)
    
    ggsave(paste0(plotPath, "figure3_", readableCt, ".png"), pPCA,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    ggsave(paste0(plotPath, "figure4_", readableCt, ".png"), pUMAP,
           width = config$figWidth,
           height = config$figHeight,
           dpi = config$figDPI,
           units = "in")
    
    message("Saving Seurat file for ", readableCt)
    saveRDS(ctObj, file = paste0(rdsPath, ct, "_seurat.rds"))
  }
}

cat("\014")
cat("\n")

message("Done!")
