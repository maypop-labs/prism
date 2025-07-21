# =============================================================================
# 09 - Boolean Regulation
#
# Infer Boolean regulation rules for each gene using binary expression data
# and structural information from the GRN. Output a list of Boolean rules.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(monocle3)
library(igraph)
library(dplyr)
library(tidyverse)
library(BoolNet)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path   <- paste0(config$rootPath, "results/monocle3/")
graphMlPath    <- paste0(config$rootPath, "results/graphml/")
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
edgesFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(degFile)) stop("Switch DEG RDS file not found: ", degFile)
if (!file.exists(edgesFile)) stop("Edge RDS file not found: ", edgesFile)
if (!file.exists(graphFile)) stop("Graph RDS file not found: ", graphFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading Switch DEGs from: ", degFile)
switchDEGs <- readRDS(degFile)
message("Loading edges from: ", edgesFile)
edges <- readRDS(edgesFile)
message("Loading graph from: ", graphFile)
g <- readRDS(graphFile)

# --- Prepare Binary Matrix and Pseudotime Order ---
matBin    <- assay(cds, "binary")
cellOrder <- order(colData(cds)$Pseudotime)
nCells    <- length(cellOrder)

winSize <- max(5, floor(config$winSizePercent * nCells)) +floor(0.5 * max(5, floor(config$winSizePercent * nCells)))
kStep <- winSize

# --- Filter Edge Table ---
if (config$grnMergeStronglyConnectedComponents) {
  gDf <- as_data_frame(g, what = "edges")
  colnames(gDf) <- c("TF", "Target", "corr", "regType")
  edges <- gDf %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = config$boolMaxRegulators) %>% ungroup()
} else {
edges <- edges %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = config$boolMaxRegulators) %>% ungroup()
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

  if (is.null(ioDf) || nrow(ioDf) < config$boolMinPairs) {
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

cat("\014")
cat("\n")

# --- Save Rules ---
if (config$saveResults) {
  message("Saving Boolean rules to: ", rulesFile)
  saveRDS(boolRules, file = rulesFile)
  
}

message("Done!")
