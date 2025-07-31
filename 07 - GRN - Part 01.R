# =============================================================================
# 07 - GRN - Part 01
#
# Build a gene regulatory network (GRN) using SCENIC from smoothed and filtered
# Monocle3 expression data and switch-like DEGs. Save the SCENIC object.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
ensureProjectDirectories(paths)
clearConsole()

# --- Load smoothed pseudotime trajectory ---
if (!dir.exists(ptPaths$monocle3SmoothedGeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3SmoothedGeneSwitches)
message("Loading Monocle3 object from: ", ptPaths$monocle3SmoothedGeneSwitches)
cds <- load_monocle_objects(directory_path = ptPaths$monocle3SmoothedGeneSwitches)

# -- load Switch DEGs ---
if (!file.exists(ptPaths$switchDegs)) stop("Switch DEG RDS file not found: ", ptPaths$switchDegs)
message("Loading switch DEGs from: ", ptPaths$switchDegs)
switchDEGS <- readRDS(ptPaths$switchDegs)

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

# --- Save Output ---
if (config$saveResults) {
  message("Saving SCENIC object to: ", ptPaths$grnPart01)
  saveRDS(scenicOptions, file = ptPaths$grnPart01)
}

message("Done!")
