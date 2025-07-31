# =============================================================================
# 08 - GRN - Part 02
#
# Use SCENIC regulons and pseudotime-filtered genes to create a signed GRN.
# Annotate interactions as activating or inhibiting via correlation.
# Optionally remove terminal nodes and merge SCCs. Save results.
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

# -- Load the SCENIC object from the previous script --
if (!file.exists(ptPaths$grnPart01)) stop("SCENIC RDS file not found: ", ptPaths$grnPart01)
message("Loading SCENIC object from: ", ptPaths$grnPart01)
scenicOptions <- readRDS(ptPaths$grnPart01)

# --- Set Up ---
sigGenes <- rownames(switchDEGs)
exprMat  <- assay(cds, "smoothed_expr")
exprMat  <- exprMat[rowSums(exprMat > 1) >= 10, ]
cellInfo <- data.frame(row.names = colnames(exprMat))

# --- Extract and Prune Regulons ---
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
}

# --- Build Graph and Prune ---
g <- graph_from_data_frame(scenicEdges, directed = TRUE)
if (config$grnRemoveTerminalNodes) {
  terminalNodes <- names(which(degree(g, mode = "out") == 0))
  g <- delete_vertices(g, terminalNodes)
}

# --- Merge Strongly Connected Components (SCCs) ---
if (config$grnMergeStronglyConnectedComponents) {
  comp <- components(g, mode = "strong")
  gMerged <- contract.vertices(g, comp$membership, vertex.attr.comb = "concat")
  g <- igraph::simplify(gMerged, remove.loops = TRUE, edge.attr.comb = "concat")  
}

# -- remove isolated nodes --
if (config$grnRemoveIsolatedNodes) {
  isolatedNodes <- names(which(degree(g, mode = "all") == 0))
  if (length(isolatedNodes) > 0) {
    g <- delete_vertices(g, isolatedNodes)
    message("Removed ", length(isolatedNodes), " isolated self-activating nodes")
  }  
}

graphPlot <- ggraph(g, layout = "fr") + # ‘fr’ = force‑directed
                    geom_edge_link(aes(colour = regType),
                    alpha = config$plotAlpha) +
                    geom_node_point(size = config$pointSize) +
                    ggtitle("Gene Regulatory Network") +
                    theme(plot.title = element_text(hjust = config$hjust)) +
                    theme_graph(background = "white") +
                    theme(legend.title = element_blank())

# --- Save Final Graph ---
if (config$saveResults) {

  message("Saving edge list to: ", ptPaths$grnPart02Edges)
  saveRDS(scenicEdges, file = ptPaths$grnPart02Edges)
  
  message("Saving final GRN to: ", ptPaths$grnPart02)
  saveRDS(g, file = ptPaths$grnPart02)
  
  message("Saving final GRN to: ", ptPaths$grnGraphml)
  igraph::write_graph(g, ptPaths$grnGraphml, format = "graphml")
  
  message("Saving plots")
  ggsave(paste0(paths$base$plots, "figure5_", cellType, "_", trajectory, ".png"), graphPlot,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
}

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
