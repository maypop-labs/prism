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
config <- initializeScript()
pathInfo <- initializeInteractivePaths()
paths <- pathInfo$paths
ensureProjectDirectories(paths)
clearConsole()

# --- Load and Annotate Seurat Objects ---
seuratList <- vector("list", length(config$donorIDs))
doubletStats <- data.frame(
  donorID = character(),
  totalCells = integer(),
  doubletsDetected = integer(),
  doubletRate = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(config$donorIDs)) {
  if (config$verbose) { message("Loading: ", config$donorIDs[i]) }
  
  counts <- Read10X(data.dir = paths$cellranger[i])
  obj <- CreateSeuratObject(counts = counts, project = config$donorIDs[i])
  obj$donorID <- config$donorIDs[i]
  obj$age <- config$ages[i]
  
  # --- Doublet Detection (per-donor) ---
  if (config$detectDoublets) {
    if (config$verbose) { 
      message("  Detecting doublets in ", config$donorIDs[i]) 
    }
    
    nCellsBefore <- ncol(obj)
    
    if (config$doubletRemovalMethod == "scDblFinder") {
      # Convert to SingleCellExperiment for scDblFinder
      sce <- as.SingleCellExperiment(obj)
      
      # Run scDblFinder
      sce <- scDblFinder(
        sce,
        dbr = config$doubletExpectedRate,
        verbose = FALSE
      )
      
      # Extract doublet classifications
      obj$doubletScore <- sce$scDblFinder.score
      obj$doubletClass <- sce$scDblFinder.class
      
      # Filter doublets
      obj <- subset(obj, subset = doubletClass == "singlet")
      
    } else if (config$doubletRemovalMethod == "DoubletFinder") {
      # DoubletFinder requires preprocessing
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, verbose = FALSE)
      obj <- ScaleData(obj, verbose = FALSE)
      obj <- RunPCA(obj, verbose = FALSE)
      obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
      
      # Run DoubletFinder parameter sweep
      sweepRes <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
      sweepStats <- summarizeSweep(sweepRes, GT = FALSE)
      bcmvn <- find.pK(sweepStats)
      pkOptimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      
      # Run DoubletFinder
      obj <- doubletFinder_v3(
        obj,
        PCs = 1:10,
        pN = 0.25,
        pK = pkOptimal,
        nExp = round(config$doubletExpectedRate * ncol(obj)),
        reuse.pANN = FALSE,
        sct = FALSE
      )
      
      # Filter doublets (column name depends on parameters)
      dfCol <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
      obj <- subset(obj, subset = obj@meta.data[[dfCol]] == "Singlet")
    }
    
    nCellsAfter <- ncol(obj)
    doubletStats <- rbind(doubletStats, data.frame(
      donorID = config$donorIDs[i],
      totalCells = nCellsBefore,
      doubletsDetected = nCellsBefore - nCellsAfter,
      doubletRate = (nCellsBefore - nCellsAfter) / nCellsBefore,
      stringsAsFactors = FALSE
    ))
    
    if (config$verbose) {
      message("  Removed ", nCellsBefore - nCellsAfter, " doublets (", 
              round(100 * (nCellsBefore - nCellsAfter) / nCellsBefore, 2), "%)")
    }
  }
  
  seuratList[[i]] <- obj
}

# Print doublet summary
if (config$detectDoublets && config$verbose) {
  message("\n--- Doublet Detection Summary ---")
  print(doubletStats)
  message("Overall doublet rate: ", 
          round(100 * sum(doubletStats$doubletsDetected) / sum(doubletStats$totalCells), 2), "%")
}

# --- Merge Seurat Objects ---
if (config$verbose) { message("\nMerging all samples into a single Seurat object") }
seuratMerged <- suppressWarnings(Reduce(function(x, y) merge(x, y, project = "All Ages and Cells"), seuratList))

# --- Quality Control ---
if (config$verbose) { message("Applying quality control filters") }
seuratMerged[["percent.mt"]] <- PercentageFeatureSet(seuratMerged, pattern = "^MT-")
seuratMerged <- subset(seuratMerged,subset = percent.mt < config$seuratMaxPercentMtDNA & nFeature_RNA > config$seuratMinNumberOfRNAFeatures)

# --- Preprocessing Pipeline ---
if (config$verbose) { message("Running Seurat v5 pipeline: normalization, PCA, clustering, and UMAP") }
seuratMerged <- NormalizeData(seuratMerged, normalization.method = "LogNormalize", scale.factor = config$seuratScaleFactor, verbose = config$verbose)
seuratMerged <- FindVariableFeatures(seuratMerged, selection.method = "vst", nfeatures = config$seuratNumberOfFeatures, verbose = config$verbose)
varFeatures <- VariableFeatures(seuratMerged)
seuratMerged <- ScaleData(seuratMerged, features = varFeatures, verbose = config$verbose)
seuratMerged <- RunPCA(seuratMerged, features = varFeatures, verbose = config$verbose)
seuratMerged <- FindNeighbors(seuratMerged, dims = 1:30, verbose = config$verbose)
seuratMerged <- FindClusters(seuratMerged, resolution = 2, cluster.name = "unintegrated_clusters", verbose = config$verbose)
seuratMerged <- RunUMAP(seuratMerged, dims = 1:30, reduction.name = "umap.unintegrated", verbose = config$verbose)

# --- Integration by Donor ---
if (config$verbose) { message("Performing donor-level integration") }
seuratMerged <- IntegrateLayers(
  object = seuratMerged,
  features = varFeatures,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = config$verbose
)
seuratMerged <- JoinLayers(seuratMerged, verbose = config$verbose)
seuratMerged <- FindNeighbors(seuratMerged, dims = 1:30, reduction = "integrated.cca", verbose = config$verbose)
seuratMerged <- FindClusters(seuratMerged, resolution = 1, verbose = config$verbose)
seuratMerged <- RunUMAP(seuratMerged, dims = 1:30, reduction = "integrated.cca", verbose = config$verbose)

# --- Plotting ---
pPCA <- DimPlot(seuratMerged,
                reduction = "pca",
                group.by = "age",
                pt.size = config$pointSize,
                stroke.size = config$strokeSize,
                alpha = config$plotAlpha) +
                ggtitle("PCA Using All Cells") +
                theme(plot.title = element_text(hjust = config$hjust)) +
                labs(x = "PC1") +
                labs(y = "PC2")

pUMAP <- DimPlot(seuratMerged,
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
  ggsave(paths$static$pcaAllCellsByAge, pPCA,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
  ggsave(paths$static$umapAllCellsByAge, pUMAP,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
  
  saveObject(seuratMerged, paths$static$seuratMerged, config, "merged Seurat object")
}

message("Done!")
