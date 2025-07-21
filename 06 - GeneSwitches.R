# =============================================================================
# 06 - GeneSwitches
#
# Perform logistic modeling of gene expression state transitions across
# pseudotime using the GeneSwitches package. Output includes filtered
# switch-like genes.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(monocle3)
library(GeneSwitches)
library(Seurat)
library(Matrix)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path   <- paste0(config$rootPath, "results/monocle3/")
plotPath       <- paste0(config$rootPath, "results/plots/")
rdsPath        <- paste0(config$rootPath, "results/rds/")
tsvPath        <- paste0(config$rootPath, "results/tsv/")
txtPath        <- paste0(config$rootPath, "results/txt/")
cellTypes      <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType       <- showCellTypeMenu(cellTypes)
trajNamesFile  <- readRDS(paste0(rdsPath, "retained_trajectories_", cellType, ".rds"))
cellTrajectory <- showTrajectoryMenu(trajNamesFile)
cdsPath        <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory, "_smoothed")
savePath       <- paste0(cdsPath, "_geneSwitches")
degFile        <- paste0(rdsPath, cellType, "_", cellTrajectory, "_degs.rds")
switchSaveFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_switch_degs.rds")
switchTsvFile  <- paste0(tsvPath, cellType, "_", cellTrajectory, "_switch_genes.tsv")

dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,  recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,  recursive = TRUE, showWarnings = FALSE)

# --- Load Required Data ---
if (!dir.exists(cdsPath)) stop("Smoothed Monocle3 directory not found: ", cdsPath)
if (!file.exists(degFile)) stop("DEG RDS file not found: ", degFile)
message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading DEG table from: ", degFile)
degTable <- readRDS(degFile)

cat("\014")
cat("\n")

# --- Filter Significant Genes ---
sigGenes <- subset(degTable, q_value < config$fdrLevel)$gene_id
message("Number of significant genes: ", length(sigGenes))
if (length(sigGenes) < config$minDEGs) stop(paste0("Fewer than ", config$minDEGs, " significant genes. Aborting."))
cds <- cds[sigGenes, ]

# --- Prepare Expression for GeneSwitches ---
message("Preparing data for GeneSwitches")
colData(cds)$Pseudotime <- pseudotime(cds)
assay(cds, "expdata") <- log1p(assay(cds, "smoothed_expr"))

# --- Run GeneSwitches Analysis ---
message("Binarizing expression")
cds <- binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesBinarizeCutoff)

message("Fitting logistic models")
cds <- find_switch_logistic_fastglm(cds, downsample = TRUE, show_warning = TRUE)

message("Filtering top switch genes")
switchGenes <- filter_switchgenes(cds, allgenes = TRUE, topnum = config$maxSwitchDEGs, r2cutoff = config$geneSwitchesR2Cutoff)

# --- Summarize switch genes for export ---
switchDf <- as.data.frame(switchGenes, stringsAsFactors = FALSE)
switchOut <- data.frame(
  geneId    = switchDf$geneID,
  direction = switchDf$direction,
  fdr       = switchDf$FDR,
  stringsAsFactors = FALSE
)
switchOut <- subset(switchOut, !is.na(fdr))
switchOut <- switchOut[order(switchOut$fdr, switchOut$geneId, na.last = TRUE), ]

cat("\014")
cat("\n")

# --- Save Results ---
if (config$saveResults) {
  message("Saving GeneSwitches Monocle3 object to: ", savePath)
  save_monocle_objects(cds = cds, directory_path = savePath, comment = cellType)

  message("Saving switch genes to: ", switchSaveFile)
  saveRDS(switchGenes, file = switchSaveFile)
  
  message("Saving switch gene report to: ", switchTsvFile)
  write.table(switchOut, file = switchTsvFile, sep  = "\t", quote = FALSE, row.names = FALSE)
}

message("Number of switch genes: ", nrow(switchOut))
message("Done!")
