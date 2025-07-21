# =============================================================================
# 07 - GRN - Part 01
#
# Build a gene regulatory network (GRN) using SCENIC from smoothed and filtered
# Monocle3 expression data and switch-like DEGs. Save the SCENIC object.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(SCENIC)
library(monocle3)
library(GENIE3)
library(Matrix)
library(dplyr)

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
cdsPath        <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile        <- paste0(rdsPath, cellType, "_", cellTrajectory, "_switch_degs.rds")
scenicSaveFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_01.rds")

dir.create(plotPath,      recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,       recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,       recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,       recursive = TRUE, showWarnings = FALSE)
dir.create(scenicOutPath, recursive = TRUE, showWarnings = FALSE)

# --- Load Monocle3 & DEG Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(degFile)) stop("Switch DEG RDS file not found: ", degFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading switch DEGs from: ", degFile)
switchDEGS <- readRDS(degFile)

cat("\014")
cat("\n")

# --- Filter and Normalize Expression Matrix ---
exprMat   <- assay(cds, "smoothed_expr")
keepGenes <- rowSums(exprMat > 1) >= 10
exprMat   <- exprMat[keepGenes, ]
cellInfo  <- data.frame(row.names = colnames(exprMat))

# --- Initialize SCENIC ---
data(defaultDbNames)
species  <- config$scenicSpecies
baseName <- paste0("motifAnnotations_", species)
# Ordered list of possible dataset names shipped in different SCENIC releases
candidates <- c(baseName, paste0(baseName, "_v10"), paste0(baseName, "_v9"), "motifAnnotations")
loaded <- FALSE
for (ds in candidates) {
  if (ds %in% data(package = "RcisTarget")$results[, "Item"]) {
    data(list = ds, package = "RcisTarget", envir = environment())
    assign(baseName, get(ds), envir = globalenv())
    loaded <- TRUE
    break
  }
}
if (!loaded) { stop("No motif annotation object found for '", species, "'. Check your RcisTarget installation or download the annotation manually.")  }

message("Initializing SCENIC")
scenicOptions <- initializeScenic(
  org           = config$scenicSpecies,
  dbDir         = config$rcisTargetPath,
  dbs           = config$scenicDBs,
  nCores        = config$cores
)

# --- Filter to DEGs and TFs ---
degGenes   <- rownames(switchDEGS)
allTFs     <- getDbTfs(scenicOptions)
unionGenes <- union(degGenes, allTFs)
exprMat    <- exprMat[intersect(rownames(exprMat), unionGenes), ]

# --- Run SCENIC Pipeline ---

message("Filtering genes and calculating co-expression")
filteredGenes <- geneFiltering(exprMat, scenicOptions)
exprMatLog    <- log2(exprMat[filteredGenes, ] + 1)
runCorrelation(exprMatLog, scenicOptions)

message("Running GENIE3")
runGenie3(exprMatLog, scenicOptions)

message("Inferring regulons and scoring")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, onlyPositiveCorr = config$grnOnlyPositiveCorr)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMatLog)

cat("\014")
cat("\n")

# --- Save Output ---
if (config$saveResults) {
  message("Saving SCENIC object to: ", scenicSaveFile)
  saveRDS(scenicOptions, file = scenicSaveFile)
}

message("Done!")
