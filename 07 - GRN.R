# =============================================================================
# 07 - GRN Construction
#
# Build a complete gene regulatory network (GRN) using SCENIC from smoothed 
# and filtered Monocle3 expression data and switch-like DEGs. Create signed
# network, prune, and save final results.
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
degGenes   <- rownames(switchDEGs)
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

# --- Extract and Prune Regulons ---
message("Building signed regulatory network")
sigGenes    <- rownames(switchDEGs)
regulons    <- loadInt(scenicOptions, "regulons")
scenicEdges <- purrr::map_dfr(names(regulons), function(tf) {
  targets   <- regulons[[tf]]
  if (length(targets) > 0) data.frame(TF = tf, Target = targets)
})
scenicEdges <- scenicEdges %>% filter(Target %in% sigGenes)

# --- Annotate Edge Sign via Correlation ---
exprMat        <- exprMat[intersect(rownames(exprMat), unique(c(scenicEdges$TF, scenicEdges$Target))), ]
validGenes     <- rownames(exprMat)
scenicEdges    <- scenicEdges %>% filter(TF %in% validGenes & Target %in% validGenes)
missingTFs     <- setdiff(scenicEdges$TF, validGenes)
missingTargets <- setdiff(scenicEdges$Target, validGenes)
if (length(missingTFs) | length(missingTargets)) {
  message("Dropped ", length(missingTFs), " TFs and ",
          length(missingTargets), " targets that were absent from exprMat after filtering.")
}

getCorrelation <- function(tf, tg) cor(exprMat[tf, ], exprMat[tg, ], method = "spearman")
scenicEdges$corr <- mapply(getCorrelation, scenicEdges$TF, scenicEdges$Target)
scenicEdges <- scenicEdges %>%
  filter(corr > config$grnPositiveThreshold |
           corr < -config$grnNegativeThreshold)
scenicEdges$regType <- ifelse(scenicEdges$corr >= 0, "Activation", "Inhibition")

# --- Add Self-Activating Loops to Orphan Targets ---
noActivatorTargets <- scenicEdges %>%
  group_by(Target) %>%
  summarize(hasActivator = any(regType == "Activation")) %>%
  filter(!hasActivator) %>% pull(Target)
if (length(noActivatorTargets) > 0) {
  selfEdges <- data.frame(
    TF      = noActivatorTargets,
    Target  = noActivatorTargets,
    regType = "Activation",
    corr    = 1
  )
  scenicEdges <- bind_rows(scenicEdges, selfEdges)
  message("Added self-activation loops for ", length(noActivatorTargets), " orphan targets")
}

# --- Build Graph and Prune ---
g <- graph_from_data_frame(scenicEdges, directed = TRUE)
message("Initial network: ", vcount(g), " nodes, ", ecount(g), " edges")

if (config$grnRemoveTerminalNodes) {
  terminalNodes <- names(which(degree(g, mode = "out") == 0))
  g <- delete_vertices(g, terminalNodes)
  message("Removed ", length(terminalNodes), " terminal nodes")
}

# --- Merge Strongly Connected Components (SCCs) ---
if (config$grnMergeStronglyConnectedComponents) {
  comp <- components(g, mode = "strong")
  gMerged <- contract.vertices(g, comp$membership, vertex.attr.comb = "concat")
  g <- igraph::simplify(gMerged, remove.loops = TRUE, edge.attr.comb = "first")
  message("Merged ", max(comp$membership), " strongly connected components")
}

# --- Remove isolated nodes ---
if (config$grnRemoveIsolatedNodes) {
  isolatedNodes <- names(which(degree(g, mode = "all") == 0))
  if (length(isolatedNodes) > 0) {
    g <- delete_vertices(g, isolatedNodes)
    message("Removed ", length(isolatedNodes), " isolated self-activating nodes")
  }  
}

message("Final network: ", vcount(g), " nodes, ", ecount(g), " edges")

# --- Create Network Plot ---
graphPlot <- ggraph(g, layout = "fr") + # 'fr' = force-directed
  geom_edge_link(aes(colour = regType),
                 alpha = config$plotAlpha) +
  geom_node_point(size = config$pointSize) +
  ggtitle("Gene Regulatory Network") +
  theme(plot.title = element_text(hjust = config$hjust)) +
  theme_graph(background = "white") +
  theme(legend.title = element_blank())

# --- Save Final Results ---
if (config$saveResults) {
  message("Saving SCENIC object to: ", ptPaths$grnPart01)
  saveRDS(scenicOptions, file = ptPaths$grnPart01)
  
  message("Saving edge list to: ", ptPaths$grnPart02Edges)
  saveRDS(scenicEdges, file = ptPaths$grnPart02Edges)
  
  message("Saving final GRN to: ", ptPaths$grnPart02)
  saveRDS(g, file = ptPaths$grnPart02)
  
  message("Saving final GRN to: ", ptPaths$grnGraphml)
  igraph::write_graph(g, ptPaths$grnGraphml, format = "graphml")
  
  message("Saving GRN plot to: ", ptPaths$grnPlot)
  ggsave(ptPaths$grnPlot, graphPlot,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
}

# --- Clean up SCENIC temporary folders ---
scenicFolders <- c("int", "output")
for (folder in scenicFolders) {
  if (dir.exists(folder)) {
    tryCatch({
      unlink(folder, recursive = TRUE)
      message("Successfully deleted SCENIC folder: ", folder)
    }, error = function(e) {
      warning("Failed to delete folder ", folder, ": ", e$message)
    })
  }
}

message("Done!")