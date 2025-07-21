# =============================================================================
# 08 - GRN - Part 02
#
# Use SCENIC regulons and pseudotime-filtered genes to create a signed GRN.
# Annotate interactions as activating or inhibiting via correlation.
# Optionally remove terminal nodes and merge SCCs. Save results.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(SCENIC)
library(monocle3)
library(dplyr)
library(Matrix)
library(igraph)
library(tidyverse)
library(ggraph)

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
graphMlFile    <- paste0(graphMlPath, cellType, "_", cellTrajectory, "_grn.graphml")
scenicFile     <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_01.rds")
edgesFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(degFile)) stop("Switch DEG RDS file not found: ", degFile)
if (!file.exists(scenicFile)) stop("SCENIC RDS file not found: ", scenicFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading Switch DEGs from: ", degFile)
switchDEGs <- readRDS(degFile)
message("Loading SCENIC object from: ", scenicFile)
scenicOptions <- readRDS(scenicFile)

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

graphPlot <- ggraph(g, layout = "fr") + # ‘fr’ = force‑directed
                    geom_edge_link(aes(colour = regType),
                    alpha = config$plotAlpha) +
                    geom_node_point(size = config$pointSize) +
                    ggtitle("Gene Regulatory Network") +
                    theme(plot.title = element_text(hjust = config$hjust)) +
                    theme_graph(background = "white") +
                    theme(legend.title = element_blank())

cat("\014")
cat("\n")

# --- Save Final Graph ---
if (config$saveResults) {

  message("Saving edge list to: ", edgesFile)
  saveRDS(scenicEdges, file = edgesFile)
  
  message("Saving final GRN to: ", graphFile)
  saveRDS(g, file = graphFile)
  
  message("Saving final GRN to: ", graphMlFile)
  igraph::write_graph(g, graphMlFile, format = "graphml")
  
  message("Saving plots")
  ggsave(paste0(plotPath, "figure5.png"), graphPlot,
         width = config$figWidth,
         height = config$figHeight,
         dpi = config$figDPI,
         units = "in")
}

message("Done!")
