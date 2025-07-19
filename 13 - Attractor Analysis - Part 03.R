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

# --- Source Helpers ---
source("functions.R")

# --- Options ---
options(warn = -1)
saveResults <- TRUE

# --- Parameters ---
cellType         <- "Keratinocytes"
cellTrajectory   <- "Y_447"
initialThreshold <- 0.25
finalThreshold   <- 0.5

# --- File Paths ---
cdsPath        <- paste0("E:/datasets/omics/skin/results/monocle3/monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
attractorsFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractors.rds")
attractorDfFile<- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractor_df.rds")
attractorScoresFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractor_df_scores.rds")
boolnetFile    <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_boolnet.rds")
boolRulesFile  <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
geneMapFile    <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_gene_map.rds")

# --- Output Paths ---
singleKDFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_single_targets_0.rds")
singleOEFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_single_targets_1.rds")
doubleKDFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_double_targets_KD_KD.rds")
doubleOEFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_double_targets_OE_OE.rds")
doubleMixFile<- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_double_targets_KD_OE.rds")

# --- Load Data ---
cds            <- load_monocle_objects(directory_path = cdsPath)
attractors     <- readRDS(attractorsFile)
attractorDf    <- readRDS(attractorDfFile)
attractorScores<- readRDS(attractorScoresFile)
boolnet        <- readRDS(boolnetFile)
boolRules      <- readRDS(boolRulesFile)
geneMap        <- readRDS(geneMapFile)

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
if (initialScore <= initialThreshold) stop("Initial score below threshold")

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
if (saveResults) {
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

if (saveResults) {
  saveRDS(doubleKD,  file = doubleKDFile)
  saveRDS(doubleOE,  file = doubleOEFile)
  saveRDS(doubleMix, file = doubleMixFile)
}

cat("Done!\n")
