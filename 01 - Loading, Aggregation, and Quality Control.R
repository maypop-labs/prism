# =============================================================================
# 01 - Loading, Aggregation, and Quality Control
#
# Load Seurat objects from CellRanger output, merge, normalize, and compute
# PCA and UMAP. Final output is an integrated Seurat object with quality
# control.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================


# --- Initialization ---
source("attractorManager.R")
source("booleanManager.R")
source("pathManager.R")
source("setupManager.R")
source("uiManager.R")

config   <- initializeScript()
pathInfo <- initializeInteractivePaths()
paths    <- pathInfo$paths
ensureProjectDirectories(paths)
clearConsole()

# --- Load and Annotate Seurat Objects ---
seuratList <- vector("list", length(config$nameVec))
for (i in seq_along(config$nameVec)) {
  message("=== Loading: ", config$nameVec[i], " ===")
  counts <- Read10X(data.dir = paths$cellranger[i])
  obj <- CreateSeuratObject(counts = counts, project = config$nameVec[i])
  obj$person <- config$nameVec[i]
  obj$age <- config$ageVec[i]
  seuratList[[i]] <- obj
}

# --- Merge Seurat Objects ---
message("Merging all samples into a single Seurat object")
mergedSeurat <- Reduce(function(x, y) merge(x, y, project = "All Ages and Cells"), seuratList)

# --- Quality Control ---
message("Applying quality control filters")
mergedSeurat[["percent.mt"]] <- PercentageFeatureSet(mergedSeurat, pattern = "^MT-")
mergedSeurat <- subset(mergedSeurat, subset = percent.mt < 15 & nFeature_RNA > 500)

# --- Preprocessing Pipeline ---
message("Running Seurat v5 pipeline: normalization, PCA, clustering, and UMAP")
mergedSeurat <- NormalizeData(mergedSeurat, normalization.method = "LogNormalize", scale.factor = 100000)
mergedSeurat <- FindVariableFeatures(mergedSeurat, selection.method = "vst", nfeatures = 2000)
varFeatures  <- VariableFeatures(mergedSeurat)
mergedSeurat <- ScaleData(mergedSeurat, features = varFeatures)
mergedSeurat <- RunPCA(mergedSeurat, features = varFeatures)
mergedSeurat <- FindNeighbors(mergedSeurat, dims = 1:30)
mergedSeurat <- FindClusters(mergedSeurat, resolution = 2, cluster.name = "unintegrated_clusters")
mergedSeurat <- RunUMAP(mergedSeurat, dims = 1:30, reduction.name = "umap.unintegrated")

# --- Integration by Donor ---
message("Performing donor-level integration")
mergedSeurat <- IntegrateLayers(
  object         = mergedSeurat,
  features       = varFeatures,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose        = config$verbose
)
mergedSeurat <- JoinLayers(mergedSeurat)
mergedSeurat <- FindNeighbors(mergedSeurat, dims = 1:30, reduction = "integrated.cca")
mergedSeurat <- FindClusters(mergedSeurat, resolution = 1)
mergedSeurat <- RunUMAP(mergedSeurat, dims = 1:30, reduction = "integrated.cca")

# --- Plotting ---
pPCA <- DimPlot(mergedSeurat,
                reduction = "pca",
                group.by = "age",
                pt.size = config$pointSize,
                stroke.size = config$strokeSize,
                alpha = config$plotAlpha) +
  ggtitle("PCA Using All Cells") +
  theme(plot.title = element_text(hjust = config$hjust)) +
  labs(x = "PC1") +
  labs(y = "PC2")

pUMAP <- DimPlot(mergedSeurat,
                 reduction = "umap.unintegrated",
                 group.by = "age",
                 pt.size = config$pointSize,
                 stroke.size = config$strokeSize,
                 alpha = config$plotAlpha) +
  ggtitle("UMAP Using All Cells") +
  theme(plot.title = element_text(hjust = config$hjust)) +
  labs(x = "UMAP1") +
  labs(y = "UMAP2")

# --- Save Results ---
if (config$saveResults) {
  
  message("Saving plots")
  ggsave(paste0(paths$static$pcaAllCells), pPCA,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
  ggsave(paste0(paths$static$umapAllCells), pUMAP,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
  
  message("Saving merged Seurat object to disk")
  saveRDS(mergedSeurat, file = paths$static$mergedSeurat)
}

message("Done!")
