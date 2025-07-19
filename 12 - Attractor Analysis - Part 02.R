# =============================================================================
# 12 - Attractor Analysis - Part 02 (Refactored)
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
options(Seurat.object.assay.version = "v5")
registerDoParallel(cores = 1)
saveResults <- TRUE

# --- Parameters ---
cellType       <- "Keratinocytes"
cellTrajectory <- "Y_447"

# --- File Paths ---
attractorsFile   <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractors.rds")
attractorDfFile  <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractor_df.rds")
boolnetFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_boolnet.rds")
outputFile       <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractor_df_scores.rds")

# --- Load Data ---
cds          <- load_monocle_objects(directory_path = paste0("E:/datasets/omics/skin/results/monocle3/monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches"))
degTable     <- readRDS(paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_switch_degs.rds"))
switchEdges  <- readRDS(paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds"))
gMerged      <- readRDS(paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02.rds"))
boolRules    <- readRDS(paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds"))
geneMap      <- readRDS(paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_gene_map.rds"))

boolnet      <- readRDS(boolnetFile)
attractors   <- readRDS(attractorsFile)
attractorDf  <- readRDS(attractorDfFile)

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

# --- Compute Final Score ---
attractorDfScores$AttractorScore <- with(attractorDfScores, Stability * BasinSize * AgeScore)
agingScore <- sum(attractorDfScores$AttractorScore)

# --- Output ---
print(attractorDfScores)
cat("\nOverall Aging Score =", agingScore, "\n")

if (saveResults) {
  message("Saving combined attractor_df_scores to: ", outputFile)
  saveRDS(attractorDfScores, file = outputFile)
}

message("Done!")
