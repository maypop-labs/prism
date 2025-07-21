# =============================================================================
# 12 - Attractor Analysis - Part 02
#
# Merge age scores and basin sizes with perturbation-based stability and entropy.
# Compute a composite attractor aging score. Save results.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(dplyr)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path        <- paste0(config$rootPath, "results/monocle3/")
graphMlPath         <- paste0(config$rootPath, "results/graphml/")
plotPath            <- paste0(config$rootPath, "results/plots/")
rdsPath             <- paste0(config$rootPath, "results/rds/")
tsvPath             <- paste0(config$rootPath, "results/tsv/")
txtPath             <- paste0(config$rootPath, "results/txt/")
cellTypes           <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType            <- showCellTypeMenu(cellTypes)
trajNamesFile       <- readRDS(paste0(rdsPath, "retained_trajectories_", cellType, ".rds"))
cellTrajectory      <- showTrajectoryMenu(trajNamesFile)
cdsPath             <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile             <- paste0(rdsPath, cellType, "_", cellTrajectory, "_switch_degs.rds")
edgesFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
boolnetFile         <- paste0(rdsPath, cellType, "_", cellTrajectory, "_boolnet.rds")
attractorsFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractors.rds")
attractorDfFile     <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractor_df.rds")
geneMapFile         <- paste0(rdsPath, cellType, "_", cellTrajectory, "_gene_map.rds")
attractorScoresFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractor_df_scores.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(boolnetFile)) stop("BoolNet RDS file not found: ", boolnetFile)
if (!file.exists(attractorsFile)) stop("Attractor RDS file not found: ", attractorsFile)
if (!file.exists(attractorDfFile)) stop("Attractor data frame RDS file not found: ", attractorDfFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading BoolNet from: ", boolnetFile)
boolnet <- readRDS(boolnetFile)
message("Loading attractors from: ", attractorsFile)
attractors <- readRDS(attractorsFile)
message("Loading attractor data frame from: ", attractorDfFile)
attractorDf <- readRDS(attractorDfFile)

# --- Compute Entropy and Stability ---
message("Computing entropy and stability metrics for attractors")
entropyDf <- computeAttractorEntropy(boolnet, attractors, nSamplesState = 1000, nPerturb = 1000)

# --- Merge Scores ---
attractorDfScores <- merge(
  x  = entropyDf,
  y  = attractorDf,
  by.x = "AttractorIndex",
  by.y = "Attractor"
)

cat("\014")
cat("\n")

# --- Compute Final Score ---
attractorDfScores$AttractorScore <- with(attractorDfScores, Stability * BasinSize * AgeScore)
agingScore <- sum(attractorDfScores$AttractorScore)

# --- Output ---
print(attractorDfScores)
cat("\nOverall Aging Score =", agingScore, "\n")

if (config$saveResults) {
  message("Saving combined attractor_df_scores to: ", attractorScoresFile)
  saveRDS(attractorDfScores, file = attractorScoresFile)
}

message("Done!")
