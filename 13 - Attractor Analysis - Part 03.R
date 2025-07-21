# =============================================================================
# 13 - Attractor Analysis - Part 03 (Refactored)
#
# Systematic in silico perturbation of genes to identify those that reduce
# the aging score. Includes both single- and double-gene perturbations.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(dplyr)
library(gmp)

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
if (!file.exists(attractorScoresFile)) stop("Attractor scores data RDS file not found: ", attractorScoresFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading BoolNet from: ", boolnetFile)
boolnet      <- readRDS(boolnetFile)
message("Loading attractors from: ", attractorsFile)
attractors <- readRDS(attractorsFile)
message("Loading attractor data frame from: ", attractorDfFile)
attractorDf <- readRDS(attractorDfFile)
message("Loading attractor scores from: ", attractorScoresFile)
attractorScores <- readRDS(attractorScoresFile)

cat("\014")
cat("\n")

# --- Output Paths ---
singleKDFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_single_targets_0.rds")
singleOEFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_single_targets_1.rds")
doubleKDFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_KD_KD.rds")
doubleOEFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_OE_OE.rds")
doubleMixFile<- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_KD_OE.rds")

# --- Helpers ---
computeAgeScore <- function(v, vYoung, vOld) {
  num <- sum((v - vYoung) * (vOld - vYoung))
  den <- sum((vOld - vYoung)^2)
  num / den
}

perturbNetworkScore <- function(network, genes, values, youngVec, oldVec, nStates = 10000) {
  net <- fixGenes(network, fixIndices = genes, values = values)
  ats <- getAttractors(net, method = "random", startStates = nStates, returnTable = TRUE)
  if (length(ats$attractors) == 0) return(NA_real_)

  entropyDf <- computeAttractorEntropy(net, ats)
  ageScores <- vapply(ats$attractors, function(a) {
    decoded <- decodeBigIntegerState(a$involvedStates[[1]], length(network$genes))
    names(decoded) <- network$genes
    gSet <- intersect(names(decoded), names(youngVec))
    computeAgeScore(decoded[gSet], youngVec[gSet], oldVec[gSet])
  }, numeric(1))

  basinSizes <- vapply(ats$attractors, function(a) a$basinSize, numeric(1))
  basinFrac  <- basinSizes / sum(basinSizes)
  entropyDf$AgeScore <- ageScores
  entropyDf$BasinSize <- basinFrac
  entropyDf$AttractorScore <- with(entropyDf, Stability * BasinSize * AgeScore)
  sum(entropyDf$AttractorScore, na.rm = TRUE)
}

# --- Define Young and Old Vectors ---
matBin <- assay(cds, "binary")
cellOrder <- order(colData(cds)$Pseudotime)
nCells <- length(cellOrder)
q20 <- floor(0.2 * nCells)
youngVec <- rowMeans(matBin[, cellOrder[1:q20]])
oldVec   <- rowMeans(matBin[, cellOrder[(nCells - q20 + 1):nCells]])

# --- Baseline Score ---
initialScore <- sum(attractorScores$AttractorScore)
cat("Initial Aging Score:", initialScore, "\n")
if (initialScore <= config$attInitialThreshold) stop("Initial score below threshold")

# --- Single-Gene Perturbations ---
genes <- intersect(boolnet$genes, names(youngVec))
singleKD <- data.frame(Gene = genes, AgingScore = NA_real_)
singleOE <- data.frame(Gene = genes, AgingScore = NA_real_)

for (gene in genes) {
  singleKD$AgingScore[singleKD$Gene == gene] <- perturbNetworkScore(boolnet, gene, 0, youngVec, oldVec)
  singleOE$AgingScore[singleOE$Gene == gene] <- perturbNetworkScore(boolnet, gene, 1, youngVec, oldVec)
}

# --- Format and Save Single Results ---
singleKD <- singleKD %>% mutate(UnperturbedScore = initialScore, Delta = AgingScore - initialScore) %>% arrange(AgingScore)
singleOE <- singleOE %>% mutate(UnperturbedScore = initialScore, Delta = AgingScore - initialScore) %>% arrange(AgingScore)

if (config$saveResults) {
  message("Saving single-gene perturbation files to: ", singleKDFile, " and ", singleOEFile)
  saveRDS(singleKD, file = singleKDFile)
  saveRDS(singleOE, file = singleOEFile)
}

# --- Double-Gene Perturbations ---
topKD <- head(singleKD$Gene, 5)
topOE <- head(singleOE$Gene, 5)

doubleKD <- expand.grid(Gene1 = topKD, Gene2 = topKD) %>% filter(Gene1 < Gene2)
doubleOE <- expand.grid(Gene1 = topOE, Gene2 = topOE) %>% filter(Gene1 < Gene2)
doubleMix <- expand.grid(KD_Gene = topKD, OE_Gene = topOE)

scorePairs <- function(df, vals) {
  df$AgingScore <- mapply(function(g1, g2) perturbNetworkScore(boolnet, c(g1, g2), vals, youngVec, oldVec), df[[1]], df[[2]])
  df$UnperturbedScore <- initialScore
  df$Delta <- df$AgingScore - initialScore
  df
}

doubleKD  <- scorePairs(doubleKD, c(0, 0))
doubleOE  <- scorePairs(doubleOE, c(1, 1))
doubleMix <- scorePairs(doubleMix, c(0, 1))

cat("\014")
cat("\n")


if (config$saveResults) {
  message("Saving double-gene perturbation files to: ", doubleKDFile, ", ", doubleOEFile, ", and ", doubleMixFile)
  saveRDS(doubleKD,  file = doubleKDFile)
  saveRDS(doubleOE,  file = doubleOEFile)
  saveRDS(doubleMix, file = doubleMixFile)
}

cat("Done!\n")
