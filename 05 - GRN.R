# =============================================================================
# 06 - GRN
#
# Build a complete gene regulatory network (GRN) with SCENIC and GENIE3 from
# smoothed and filtered Monocle3 expression data and switch-like genes. Create
# signed network, prune, and save final results.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/grnManager.R")
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
if (!dir.exists(ptPaths$monocle3GeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
if (config$verbose) { message("Loading Monocle3 object from: ", ptPaths$monocle3GeneSwitches) }
cds <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)

# -- load Switch Genes ---
if (!file.exists(ptPaths$geneSwitches)) stop("Switch genes RDS file not found: ", ptPaths$geneSwitches)
if (config$verbose) { message("Loading switch genes from: ", ptPaths$geneSwitches) }
switchGenes <- readRDS(ptPaths$geneSwitches)


# --- Filter and Normalize Expression Matrix ---
switchGeneNames <- rownames(switchGenes)
exprMat         <- assay(cds, "counts")
keepGenes       <- rowSums(exprMat > 1) >= 10
exprMat         <- exprMat[keepGenes, ]
cellInfo        <- data.frame(row.names = colnames(exprMat))

# --- Initialize SCENIC ---
data(defaultDbNames)
species  <- config$grnScenicSpecies
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

if (config$verbose) { message("Initializing SCENIC") }
scenicOptions <- initializeScenic(
  org    = config$grnScenicSpecies,
  dbDir  = config$rcisTargetPath,
  dbs    = config$grnScenicDBs,
  nCores = config$cores
)

# --- Filter to switch genes and TFs ---
allTFs     <- getDbTfs(scenicOptions)
unionGenes <- union(switchGeneNames, allTFs)
exprMat    <- exprMat[intersect(rownames(exprMat), unionGenes), ]

# --- Run SCENIC Pipeline ---
if (config$verbose) { message("Filtering genes and calculating co-expression") }
filteredGenes <- geneFiltering(exprMat, scenicOptions)
exprMatLog    <- log2(exprMat[filteredGenes, ] + 1)
runCorrelation(exprMatLog, scenicOptions)

if (config$verbose) { message("Running GENIE3") }
runGenie3(exprMatLog, scenicOptions)

if (config$verbose) { message("Inferring regulons and scoring") }
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, onlyPositiveCorr = config$grnOnlyPositiveCorr)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMatLog)

# --- Extract and Prune Regulons ---
if (config$verbose) { message("Building signed regulatory network") }
regulons    <- loadInt(scenicOptions, "regulons")
scenicEdges <- purrr::map_dfr(names(regulons), function(tf) {
  targets   <- regulons[[tf]]
  if (length(targets) > 0) data.frame(TF = tf, Target = targets)
})
scenicEdges <- scenicEdges %>% filter(Target %in% switchGeneNames)

# --- Annotate Edge Sign via Correlation ---
exprMat        <- exprMat[intersect(rownames(exprMat), unique(c(scenicEdges$TF, scenicEdges$Target))), ]
validGenes     <- rownames(exprMat)
scenicEdges    <- scenicEdges %>% filter(TF %in% validGenes & Target %in% validGenes)
missingTFs     <- setdiff(scenicEdges$TF, validGenes)
missingTargets <- setdiff(scenicEdges$Target, validGenes)
if (length(missingTFs) | length(missingTargets)) {
  if (config$verbose) { message("Dropped ", length(missingTFs), " TFs and ", length(missingTargets), " targets that were absent from exprMat after filtering.") }
}

getCorrelation      <- function(tf, tg) cor(exprMat[tf, ], exprMat[tg, ], method = "spearman")
scenicEdges$corr    <- mapply(getCorrelation, scenicEdges$TF, scenicEdges$Target)
scenicEdges         <- scenicEdges %>% filter(corr > config$grnPositiveThreshold | corr < -config$grnNegativeThreshold)
scenicEdges$regType <- ifelse(scenicEdges$corr >= 0, "Activation", "Inhibition")

# --- Build Graph and Prune ---
g <- graph_from_data_frame(scenicEdges, directed = TRUE)

# --- Remove isolated nodes ---
if (config$grnRemoveIsolatedNodes) {
  isolatedNodes <- names(which(degree(g, mode = "all") == 0))
  if (length(isolatedNodes) > 0) {
    g <- delete_vertices(g, isolatedNodes)
    if (config$verbose) { message("Removed ", length(isolatedNodes), " isolated nodes") }
  }  
}

if (config$verbose) { message("Final network: ", vcount(g), " nodes, ", ecount(g), " edges") }

# --- Save results ---
if (config$saveResults) {
  saveGrnOutputSet(
    graph = g,
    rdsPath = ptPaths$grn,
    plotPath = ptPaths$grnPlot,
    graphmlPath = ptPaths$grnGraphml,
    title = "Gene Regulatory Network",
    config = config,
    verbose = config$verbose
  )

  if (config$verbose) { message("Saving GRN RDS file to: ", ptPaths$grn) }
  saveRDS(g, file = ptPaths$grn)
  
  if (config$verbose) { message("Saving edge list RDS file to: ", ptPaths$grnEdges) }
  saveRDS(scenicEdges, file = ptPaths$grnEdges)

}

# --- Clean up SCENIC temporary folders ---
scenicFolders <- c("int", "output")
for (folder in scenicFolders) {
  if (dir.exists(folder)) {
    tryCatch({
      unlink(folder, recursive = TRUE)
      if (config$verbose) { message("Successfully deleted SCENIC folder: ", folder) }
    }, error = function(e) {
      if (config$verbose) { warning("Failed to delete folder ", folder, ": ", e$message) }
    })
  }
}

message("Done!")
