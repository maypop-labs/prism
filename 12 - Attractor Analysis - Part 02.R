# =============================================================================
# 12 - Attractor Analysis - Part 02
#
# Merge age scores and basin sizes with perturbation-based stability and entropy.
# Compute a composite attractor aging score. Save results.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(dplyr)
library(gmp)
library(progressr)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
cl <- makePSOCKcluster(config$cores)
clusterExport(cl, c("decodeBigIntegerState", "shannonEntropy"))
registerDoParallel(cl)

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

cat("\014")

# =============================================================================
# Compute attractor entropy and stability
# =============================================================================

message("Computing entropy (", config$nSamplesState, " states per attractor, ", config$nPerturb, " perturbations each)")

handlers("txtprogressbar")

entropyDf <- switch(config$entropyMethod,
                    shuffle = entropyDf <- computeAttractorEntropy(
                      boolnet       = boolnet,
                      attractors    = attractors,
                      nSamplesState = config$nSamplesState,
                      nPerturb      = config$nPerturb,
                      showProgress  = TRUE),
                    bitflip = entropyDf <- computeAttractorEntropy_bitflip(
                      boolnet       = boolnet,
                      attractors    = attractors,
                      nSamplesState = config$nSamplesState,
                      nPerturb      = config$nPerturb,
                      showProgress  = TRUE)
)

# =============================================================================
# Merge with basin size and age scores
# =============================================================================

message("Merging scores…")
attractorDfScores <- merge(
  x  = entropyDf,
  y  = attractorDf,
  by.x = "AttractorIndex",
  by.y = "Attractor",
  all.x = TRUE)

# =============================================================================
# Composite aging score
# =============================================================================

attractorDfScores$AttractorScore <- with(attractorDfScores,
                                         Stability * BasinSize * AgeScore)
agingScore <- sum(attractorDfScores$AttractorScore, na.rm = TRUE)

# =============================================================================
# Output
# =============================================================================

print(attractorDfScores)
cat("\nOverall Aging Score =", agingScore, "\n")

if (isTRUE(config$saveResults)) {
  message("Saving combined attractor_df_scores to: ", attractorScoresFile)
  saveRDS(attractorDfScores, file = attractorScoresFile)
}

message("All done!")
