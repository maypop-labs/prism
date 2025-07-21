# =============================================================================
# 10 - BoolNet
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
rulesFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
boolnetFile    <- paste0(rdsPath, cellType, "_", cellTrajectory, "_boolnet.rds")
attractorsFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractors.rds")
geneMapFile    <- paste0(rdsPath, cellType, "_", cellTrajectory, "_gene_map.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(rulesFile)) stop("Boolean rules RDS file not found: ", graphFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading Boolean Rules from: ", rulesFile)
boolRules <- readRDS(rulesFile)

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
  startStates = config$boolNetStartStates,
  returnTable = TRUE,
  type = "synchronous"
)

cat("\014")
cat("\n")

# --- Save Outputs ---
if (config$saveResults) {
  message("Saving BoolNet object to: ", boolnetFile)
  saveRDS(boolnet, file = boolnetFile)

  message("Saving attractors to: ", attractorsFile)
  saveRDS(attractors, file = attractorsFile)

  message("Saving gene mapping to: ", geneMapFile)
  geneMap <- generateSanitizedGeneMapping(rownames(cds))
  saveRDS(geneMap, file = geneMapFile)
}

message("Done!")
