# =============================================================================
# 10 - BoolNet (Refactored)
#
# Convert Boolean rules into a BoolNet object, compute attractors,
# and save the BoolNet network, attractors, and sanitized gene name mapping.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(igraph)
library(dplyr)

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

cdsPath        <- paste0("E:/datasets/omics/skin/results/monocle3/monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile        <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_switch_degs.rds")
edgesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile      <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_Boolean_Rules.rds")

boolnetSaveFile    <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_boolnet.rds")
attractorsSaveFile <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_attractors.rds")
geneMapSaveFile    <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_gene_map.rds")

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found")
if (!file.exists(degFile)) stop("DEG file not found")
if (!file.exists(edgesFile)) stop("Edge file not found")
if (!file.exists(graphFile)) stop("Graph file not found")
if (!file.exists(rulesFile)) stop("Boolean rule file not found")

cds        <- load_monocle_objects(directory_path = cdsPath)
degTable   <- readRDS(degFile)
scenicEdges <- readRDS(edgesFile)
gMerged    <- readRDS(graphFile)
boolRules  <- readRDS(rulesFile)

# --- Build BoolNet Rule Table ---
message("Building BoolNet rule table")
ruleLines <- c("targets,factors")
for (gene in names(boolRules)) {
  safeGene <- sanitizeGeneName(gene)
  safeRegs <- sapply(boolRules[[gene]]$regulators, sanitizeGeneName, USE.NAMES = FALSE)
  pattern  <- boolRules[[gene]]$outPattern
  ruleStr  <- makeBoolNetRule(safeGene, safeRegs, pattern)
  ruleLines <- c(ruleLines, ruleStr)
}

# --- Convert to BoolNet Object ---
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)
boolnet <- loadNetwork(tempFile)
boolnet$type <- "synchronous"

# --- Compute Attractors ---
message("Computing attractors with BoolNet")
attractors <- getAttractors(
  boolnet,
  method = "random",
  startStates = 100000,
  returnTable = TRUE,
  type = "synchronous"
)

# --- Save Outputs ---
if (saveResults) {
  message("Saving BoolNet object")
  saveRDS(boolnet, file = boolnetSaveFile)

  message("Saving attractors")
  saveRDS(attractors, file = attractorsSaveFile)

  message("Saving gene mapping")
  geneMap <- generateSanitizedGeneMapping(rownames(cds))
  saveRDS(geneMap, file = geneMapSaveFile)
}

message("Done!")
