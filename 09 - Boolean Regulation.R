# =============================================================================
# 09 - Boolean Regulation (Refactored)
#
# Infer Boolean regulation rules for each gene using binary expression data
# and structural information from the GRN. Output a list of Boolean rules.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(monocle3)
library(igraph)
library(dplyr)
library(tidyverse)
library(BoolNet)

# --- Source Helpers ---
source("functions.R")

# --- Options ---
options(warn = -1)
options(Seurat.object.assay.version = "v5")
registerDoParallel(cores = 1)
saveResults <- TRUE
useMergedGraph <- FALSE

# --- Parameters ---
cellType       <- "Keratinocytes"
cellTrajectory <- "Y_447"
maxRegulators  <- 3
minPairs       <- 5

cdsPath        <- paste0("E:/datasets/omics/skin/results/monocle3/monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile        <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_switch_degs.rds")
graphFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
edgesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
ruleSaveFile   <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds")

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("CDS directory not found")
if (!file.exists(degFile)) stop("DEG RDS not found")
if (!file.exists(edgesFile)) stop("Edge RDS not found")
if (!file.exists(graphFile)) stop("Graph RDS not found")

message("Loading CDS")
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading DEGs")
degTable <- readRDS(degFile)
message("Loading edges")
edges <- readRDS(edgesFile)
message("Loading graph")
gMerged <- readRDS(graphFile)

# --- Prepare Binary Matrix and Pseudotime Order ---
matBin <- assay(cds, "binary")
cellOrder <- order(colData(cds)$Pseudotime)
nCells <- length(cellOrder)

winSize <- max(5, floor(0.10 * nCells)) + floor(0.5 * max(5, floor(0.10 * nCells)))
kStep <- winSize

# --- Filter Edge Table ---
if (useMergedGraph) {
  mergedDf <- as_data_frame(gMerged, what = "edges")
  colnames(mergedDf) <- c("TF", "Target", "corr", "regType")
  edges <- mergedDf %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = maxRegulators) %>% ungroup()
} else {
  edges <- edges %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = maxRegulators) %>% ungroup()
}

# --- Filter to Genes Present in Expression Matrix ---
geneList <- rownames(matBin)
adjTable <- edges %>% filter(TF %in% geneList & Target %in% geneList) %>% select(TF, Target)

# --- Infer Boolean Rules ---
message("Inferring Boolean rules for ", length(unique(adjTable$Target)), " target genes")
boolRules <- list()

geneIndex <- 1
for (gene in unique(adjTable$Target)) {
  message("[", geneIndex, "/", length(unique(adjTable$Target)), "] Gene: ", gene)
  geneIndex <- geneIndex + 1

  regulators <- unique(adjTable$TF[adjTable$Target == gene])
  if (length(regulators) == 0) {
    message("No regulators. Skipping.")
    next
  }

  ioDf <- makeInputOutputPairs(
    targetGene = gene,
    regulators = regulators,
    matBin     = matBin,
    cellOrder  = cellOrder,
    k          = kStep
  )

  if (is.null(ioDf) || nrow(ioDf) < minPairs) {
    message("Too few input-output pairs. Skipping.")
    next
  }

  res <- findBestBooleanRules(ioDf)
  pattern <- combineBooleanFunctionsByOr(res$bestFns)

  boolRules[[gene]] <- list(
    regulators  = regulators,
    bestScore   = res$score,
    outPattern  = pattern
  )
}

# --- Save Rules ---
if (saveResults) {
  message("Saving Boolean rules to: ", ruleSaveFile)
  saveRDS(boolRules, file = ruleSaveFile)
}

message("Done!")
