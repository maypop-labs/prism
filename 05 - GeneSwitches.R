# =============================================================================
# 05 - GeneSwitches
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
if (config$verbose) { message("Loading Monocle3 object from: ", ptPaths$monocle3Smoothed) }
cds <- load_monocle_objects(directory_path = ptPaths$monocle3Smoothed)

# --- Prepare Expression for GeneSwitches ---
if (config$verbose) { message("Preparing data for GeneSwitches") }
colData(cds)$Pseudotime <- pseudotime(cds)
assay(cds, "expdata") <- log1p(assay(cds, "smoothed_expr"))

# --- Run GeneSwitches Analysis ---
if (config$verbose) { message("Binarizing expression") }
cds <- binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesBinarizeCutoff)

# Filter out genes with zero variance after binarization
binary_assay <- assay(cds, "binary")
gene_variance <- apply(binary_assay, 1, var, na.rm = TRUE)
variable_genes <- rownames(cds)[gene_variance > 0]
cds <- cds[variable_genes, ]
if (config$verbose) { message("Retained ", length(variable_genes), " genes with binary variation") }

if (config$verbose) { message("Fitting logistic models") }
cds <- find_switch_logistic_fastglm(cds, downsample = TRUE, show_warning = TRUE)

if (config$verbose) { message("Filtering top switch genes") }
switchGenes <- filter_switchgenes(cds, allgenes = TRUE, topnum = config$geneSwitchesMaxGenes, r2cutoff = config$geneSwitchesR2Cutoff)

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
  if (config$verbose) { message("Saving GeneSwitches Monocle3 object to: ", ptPaths$monocle3GeneSwitches) }
  save_monocle_objects(cds = cds, directory_path =  ptPaths$monocle3GeneSwitches, comment = cellType)

  if (config$verbose) { message("Saving switch genes to: ", ptPaths$geneSwitches) }
  saveRDS(switchGenes, file = ptPaths$geneSwitches)
  
  if (config$verbose) { message("Saving switch gene report to: ", ptPaths$geneSwitchesTsv) }
  write.table(switchOut, file = ptPaths$geneSwitchesTsv, sep  = "\t", quote = FALSE, row.names = FALSE)
}

if (config$verbose) { message("Number of switch genes: ", nrow(switchOut)) }
message("Done!")
