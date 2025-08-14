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

# --- Load and Annotate Seurat Objects ---
seuratList <- vector("list", length(config$donorIDs))
for (i in seq_along(config$donorIDs)) {
  if (config$verbose) { message("Loading: ", config$donorIDs[i]) }
  counts          <- Read10X(data.dir = paths$cellranger[i])
  obj             <- CreateSeuratObject(counts = counts, project = config$donorIDs[i])
  obj$donorID     <- config$donorIDs[i]
  obj$age         <- config$ages[i]
  seuratList[[i]] <- obj
}

# --- Merge Seurat Objects ---
if (config$verbose) { message("Merging all samples into a single Seurat object") }
mergedSeurat <- Reduce(function(x, y) merge(x, y, project = "All Ages and Cells"), seuratList)

# --- Quality Control ---
if (config$verbose) { message("Applying quality control filters") }
mergedSeurat[["percent.mt"]] <- PercentageFeatureSet(mergedSeurat, pattern = "^MT-")
mergedSeurat <- subset(mergedSeurat,subset = percent.mt < config$seuratMaxPercentMtDNA & nFeature_RNA > config$seuratMinNumberOfRNAFeatures)

# --- Preprocessing Pipeline ---
if (config$verbose) { message("Running Seurat v5 pipeline: normalization, PCA, clustering, and UMAP") }
mergedSeurat <- NormalizeData(mergedSeurat, normalization.method = "LogNormalize", scale.factor = config$seuratScaleFactor)
mergedSeurat <- FindVariableFeatures(mergedSeurat, selection.method = "vst", nfeatures = config$seuratNumberOfFeatures)
varFeatures  <- VariableFeatures(mergedSeurat)
mergedSeurat <- ScaleData(mergedSeurat, features = varFeatures)
mergedSeurat <- RunPCA(mergedSeurat, features = varFeatures)
mergedSeurat <- FindNeighbors(mergedSeurat, dims = 1:30)
mergedSeurat <- FindClusters(mergedSeurat, resolution = 2, cluster.name = "unintegrated_clusters")
mergedSeurat <- RunUMAP(mergedSeurat, dims = 1:30, reduction.name = "umap.unintegrated")

# --- Integration by Donor ---
if (config$verbose) { message("Performing donor-level integration") }
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
  
  if (config$verbose) { message("Saving plots") }
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
  
  if (config$verbose) { message("Saving merged Seurat object to disk") }
  saveRDS(mergedSeurat, file = paths$static$mergedSeurat)
}

message("Done!")
