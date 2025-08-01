# =============================================================================
# 08 - Boolean Regulation
#
# Infer Boolean regulation rules for each gene using binary expression data
# and structural information from the GRN. Output a list of Boolean rules.
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
switchDEGs <- readRDS(ptPaths$switchDegs)

# -- Load GRN edge list ---
if (!file.exists(ptPaths$grnPart02Edges)) stop("Edges RDS file not found: ", ptPaths$grnPart02Edges)
message("Loading edges from: ", ptPaths$grnPart02Edges)
edges <- readRDS(ptPaths$grnPart02Edges)

if (!file.exists(ptPaths$grnPart02)) stop("Final GRN RDS file not found: ", ptPaths$grnPart02)
message("Loading final GRN from: ", ptPaths$grnPart02)
g <- readRDS(ptPaths$grnPart02)
stop()
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

# --- Save Rules ---
if (config$saveResults) {
  message("Saving Boolean rules to: ", ptPaths$booleanRules)
  saveRDS(boolRules, file = ptPaths$booleanRules)
  
}

message("Done!")
