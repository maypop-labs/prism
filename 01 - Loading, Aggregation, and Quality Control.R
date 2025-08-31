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
source("managers/pathManager.R")
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
seuratMerged <- suppressWarnings(Reduce(function(x, y) merge(x, y, project = "All Ages and Cells"), seuratList))

# --- Quality Control ---
if (config$verbose) { message("Applying quality control filters") }
seuratMerged[["percent.mt"]] <- PercentageFeatureSet(seuratMerged, pattern = "^MT-")
seuratMerged <- subset(seuratMerged,subset = percent.mt < config$seuratMaxPercentMtDNA & nFeature_RNA > config$seuratMinNumberOfRNAFeatures)

# --- Preprocessing Pipeline ---
if (config$verbose) { message("Running Seurat v5 pipeline: normalization, PCA, clustering, and UMAP") }
seuratMerged <- NormalizeData(seuratMerged, normalization.method = "LogNormalize", scale.factor = config$seuratScaleFactor, verbose = config$verbose)
seuratMerged <- FindVariableFeatures(seuratMerged, selection.method = "vst", nfeatures = config$seuratNumberOfFeatures, verbose = config$verbose)
varFeatures  <- VariableFeatures(seuratMerged)
seuratMerged <- ScaleData(seuratMerged, features = varFeatures, verbose = config$verbose)
seuratMerged <- RunPCA(seuratMerged, features = varFeatures, verbose = config$verbose)
seuratMerged <- FindNeighbors(seuratMerged, dims = 1:30, verbose = config$verbose)
seuratMerged <- FindClusters(seuratMerged, resolution = 2, cluster.name = "unintegrated_clusters", verbose = config$verbose)
seuratMerged <- RunUMAP(seuratMerged, dims = 1:30, reduction.name = "umap.unintegrated", verbose = config$verbose)

# --- Integration by Donor ---
if (config$verbose) { message("Performing donor-level integration") }
seuratMerged <- IntegrateLayers(
  object         = seuratMerged,
  features       = varFeatures,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose        = config$verbose
)
seuratMerged <- JoinLayers(seuratMerged, verbose = config$verbose)
seuratMerged <- FindNeighbors(seuratMerged, dims = 1:30, reduction = "integrated.cca", verbose = config$verbose)
seuratMerged <- FindClusters(seuratMerged, resolution = 1, verbose = config$verbose)
seuratMerged <- RunUMAP(seuratMerged, dims = 1:30, reduction = "integrated.cca", verbose = config$verbose)

# --- Plotting ---
pPCA <- DimPlot(seuratMerged,
                reduction   = "pca",
                group.by    = "age",
                pt.size     = config$pointSize,
                stroke.size = config$strokeSize,
                alpha       = config$plotAlpha) +
                ggtitle("PCA Using All Cells") +
                theme(plot.title = element_text(hjust = config$hjust)) +
                labs(x = "PC1") +
                labs(y = "PC2")

pUMAP <- DimPlot(seuratMerged,
                 reduction   = "umap.unintegrated",
                 group.by    = "age",
                 pt.size     = config$pointSize,
                 stroke.size = config$strokeSize,
                 alpha       = config$plotAlpha) +
                 ggtitle("UMAP Using All Cells") +
                 theme(plot.title = element_text(hjust = config$hjust)) +
                 labs(x = "UMAP1") +
                 labs(y = "UMAP2")

# --- Save Results ---
if (config$saveResults) {
  
  if (config$verbose) { message("Saving plots") }
  ggsave(paths$static$pcaAllCellsByAge, pPCA,
         width  = config$figWidth,
         height = config$figHeight,
         dpi    = config$figDPI,
         units  = "in")
  ggsave(paths$static$umapAllCellsByAge, pUMAP,
         width  = config$figWidth,
         height = config$figHeight,
         dpi    = config$figDPI,
         units  = "in")
  
  saveMergedSeurat(seuratMerged, paths, config)
}

message("Done!")
