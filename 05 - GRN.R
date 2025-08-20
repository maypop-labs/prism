# =============================================================================
# 05 - GRN
#
# Build a complete gene regulatory network (GRN) with integration of
# SCENIC/GENIE3 biological priors. Uses composite ranking, motif-aware signing,
# and metadata extraction.
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

# --- Load pseudotime trajectory ---
if (!dir.exists(ptPaths$monocle3GeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
if (config$verbose) { message("Loading Monocle3 object from: ", ptPaths$monocle3GeneSwitches) }
cds <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)

# --- Load Switch Genes ---
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

# Convert sparse matrix to dense matrix for SCENIC compatibility
if (methods::is(exprMat, "dgCMatrix") || methods::is(exprMat, "Matrix")) {
  if (config$verbose) { message("Converting sparse matrix to dense matrix for SCENIC") }
  exprMat <- as.matrix(exprMat)
}

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
if (config$verbose) { message("Building enhanced regulatory network with biological priors") }
regulons    <- loadInt(scenicOptions, "regulons")
scenicEdges <- purrr::map_dfr(names(regulons), function(tf) {
  targets   <- regulons[[tf]]
  if (length(targets) > 0) data.frame(TF = tf, Target = targets)
})
scenicEdges <- scenicEdges %>% filter(Target %in% switchGeneNames)

# --- Extract SCENIC Metadata ---
scenicEdges <- extractScenicMetadata(scenicOptions, scenicEdges, config$verbose)

# --- Annotate Edge Sign via Correlation and Filter ---
exprMat        <- exprMat[intersect(rownames(exprMat), unique(c(scenicEdges$TF, scenicEdges$Target))), ]
validGenes     <- rownames(exprMat)
scenicEdges    <- scenicEdges %>% filter(TF %in% validGenes & Target %in% validGenes)
missingTFs     <- setdiff(scenicEdges$TF, validGenes)
missingTargets <- setdiff(scenicEdges$Target, validGenes)
if (length(missingTFs) | length(missingTargets)) {
  if (config$verbose) { message("Dropped ", length(missingTFs), " TFs and ", length(missingTargets), " targets that were absent from exprMat after filtering.") }
}
getCorrelation   <- function(tf, tg) cor(exprMat[tf, ], exprMat[tg, ], method = "spearman")
scenicEdges$corr <- mapply(getCorrelation, scenicEdges$TF, scenicEdges$Target)

# --- Filter edges ---
scenicEdges <- assignMotifAwareSigns(scenicEdges, config, config$verbose)
scenicEdges <- computeCompositeRanking(scenicEdges, config, config$verbose)
keepEdges   <- with(scenicEdges, corr > config$grnPositiveThreshold | corr < -config$grnNegativeThreshold | priorStrength > config$grnPriorThreshold)
scenicEdges <- scenicEdges[keepEdges, ]

if (config$verbose) {
  message("Retained ", nrow(scenicEdges), " edges after filtering")
  corrOnlyKeep  <- sum(with(scenicEdges, abs(corr) > config$grnPositiveThreshold & priorStrength <= config$grnPriorThreshold))
  priorOnlyKeep <- sum(with(scenicEdges, abs(corr) <= config$grnPositiveThreshold & priorStrength > config$grnPriorThreshold))
  bothKeep      <- sum(with(scenicEdges, abs(corr) > config$grnPositiveThreshold & priorStrength > config$grnPriorThreshold))
  message("  - Kept by correlation only: ", corrOnlyKeep)
  message("  - Kept by prior only: ", priorOnlyKeep) 
  message("  - Kept by both: ", bothKeep)
}

# --- Build Graph ---
g <- graph_from_data_frame(scenicEdges, directed = TRUE)
if (config$verbose) { message("Initial network: ", vcount(g), " nodes, ", ecount(g), " edges") }

# --- Clean up ---
if (config$grnRemoveIsolatedNodes) {
  isolatedNodes <- names(which(degree(g, mode = "all") == 0))
  if (length(isolatedNodes) > 0) {
    g <- delete_vertices(g, isolatedNodes)
    # Update edge list
    scenicEdges <- scenicEdges %>% 
      filter(TF %in% igraph::V(g)$name & Target %in% igraph::V(g)$name)
    if (config$verbose) { message("Removed ", length(isolatedNodes), " isolated nodes") }
  }  
}

if (config$verbose) { message("Final network: ", vcount(g), " nodes, ", ecount(g), " edges") }

# --- Save Final Enhanced Results ---
if (config$saveResults) {
  saveEnhancedGrnOutputSet(
    graph       = g,
    edges       = scenicEdges,
    rdsPath     = ptPaths$grn,
    plotPath    = ptPaths$grnPlot,
    graphmlPath = ptPaths$grnGraphml,
    edgesPath   = ptPaths$grnEdges,
    title       = "Gene Regulatory Network (GRN)",
    config      = config,
    verbose     = config$verbose
  )
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