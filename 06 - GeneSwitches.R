# =============================================================================
# 06 - GeneSwitches
#
# Perform logistic modeling of gene expression state transitions across
# pseudotime using the GeneSwitches package. Output includes filtered
# switch-like genes.
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
if (!dir.exists(ptPaths$monocle3Smoothed)) stop("Smoothed Monocle3 directory not found: ", ptPaths$monocle3Smoothed)
message("Loading Monocle3 object from: ", ptPaths$monocle3Smoothed)
cds <- load_monocle_objects(directory_path = ptPaths$monocle3Smoothed)

# --- Load DEG file ---
if (!file.exists(ptPaths$degs)) stop("DEG RDS file not found: ", ptPaths$degs)
message("Loading DEG table from: ", ptPaths$degs)
degTable <- readRDS(ptPaths$degs)

# --- Filter Significant Genes ---
sigGenes <- subset(degTable, q_value < config$fdrLevel)$gene_id
message("Number of significant genes: ", length(sigGenes))
if (length(sigGenes) < config$minDEGs) stop(paste0("Fewer than ", config$minDEGs, " significant genes. Aborting."))
cds <- cds[sigGenes, ]

# --- Prepare Expression for GeneSwitches ---
message("Preparing data for GeneSwitches")
colData(cds)$Pseudotime <- pseudotime(cds)
assay(cds, "expdata") <- log1p(assay(cds, "smoothed_expr"))

# --- Run GeneSwitches Analysis ---
message("Binarizing expression")
cds <- binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesBinarizeCutoff)

message("Fitting logistic models")
cds <- find_switch_logistic_fastglm(cds, downsample = TRUE, show_warning = TRUE)

message("Filtering top switch genes")
switchGenes <- filter_switchgenes(cds, allgenes = TRUE, topnum = config$maxSwitchDEGs, r2cutoff = config$geneSwitchesR2Cutoff)

# --- Summarize switch genes for export ---
switchDf <- as.data.frame(switchGenes, stringsAsFactors = FALSE)
switchOut <- data.frame(
  geneId    = switchDf$geneID,
  direction = switchDf$direction,
  fdr       = switchDf$FDR,
  stringsAsFactors = FALSE
)
switchOut <- subset(switchOut, !is.na(fdr))
switchOut <- switchOut[order(switchOut$fdr, switchOut$geneId, na.last = TRUE), ]

# --- Save Results ---
if (config$saveResults) {
  message("Saving GeneSwitches Monocle3 object to: ", ptPaths$monocle3SmoothedGeneSwitches)
  save_monocle_objects(cds = cds, directory_path =  ptPaths$monocle3SmoothedGeneSwitches, comment = cellType)

  message("Saving switch genes to: ", ptPaths$switchDegs)
  saveRDS(switchGenes, file = ptPaths$switchDegs)
  
  message("Saving switch gene report to: ", ptPaths$degsTsv)
  write.table(switchOut, file = ptPaths$degsTsv, sep  = "\t", quote = FALSE, row.names = FALSE)
}

message("Number of switch genes: ", nrow(switchOut))
message("Done!")
