# =============================================================================
# 04 - GeneSwitches
#
# Perform logistic modeling of gene expression state transitions across
# pseudotime using the GeneSwitches package on raw data.
# Output includes filtered switch-like genes.
# =============================================================================

# --- Initialization ---
source("managers/pathManager.R")
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

# --- Load Data ---
cds <- loadMonocle3(ptPaths$monocle3, config, "pseudotime trajectory")

# --- Prepare Data for GeneSwitches ---
if (config$verbose) { message("Preparing data for GeneSwitches") }
colData(cds)$Pseudotime <- pseudotime(cds)
assay(cds, "expdata")   <- log1p(assay(cds, "counts"))
expdata                 <- assay(cds, "expdata")
pt                      <- colData(cds)$Pseudotime
detect_frac             <- rowMeans(expdata > 0, na.rm = TRUE)
keep                    <- (detect_frac >= config$geneSwitchesMinDetect)

if (config$verbose) message("Prefiltering by detect ≥ ", config$geneSwitchesMinDetect)
if (config$verbose) message("Genes kept pre-binarization: ", sum(keep), " / ", nrow(expdata))
if (sum(keep) == 0L) { stop("No genes passed the prefilter.") }
cds     <- cds[keep, ]
expdata <- assay(cds, "expdata")

# =============================================================================
# --- Core GeneSwitches Pipeline ---
if (config$geneSwitchesFixedCutoff) {

  # Fixed cutoff path
  if (config$verbose) { message("Running GeneSwitches Pipeline (fixed cutoff = ", config$geneSwitchesCutOff, ")") }
  cds <- binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesCutOff, ncores = config$cores)

} else {
  
  # Mixture model path (with fallback)
  if (config$verbose) message("Running GeneSwitches Pipeline")
  # Make sure expdata is a numeric matrix
  exp <- assay(cds, "expdata")
  if (!is.matrix(exp)) exp <- as.matrix(exp)
  storage.mode(exp) <- "double"
  exp[!is.finite(exp)] <- 0
  assay(cds, "expdata") <- exp
  cds <- tryCatch(
    { binarize_exp(cds, fix_cutoff = FALSE, ncores = config$cores) },
    error = function(e) {
      warning("Mixture binarization failed: ", conditionMessage(e), " — falling back to fixed cutoff.")
      binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesCutOff, ncores = config$cores)
    }
  )
}

# Fit logistic models (avoid downsampling in sparse data)
cds <- find_switch_logistic_fastglm(cds, downsample = FALSE, show_warnings = TRUE)

# Filter/collect switch genes (lenient; all genes considered)
switchGenes <- filter_switchgenes(cds, allgenes = TRUE)

# =============================================================================

# --- Summarize Switch Genes for Export ---
switchDf  <- as.data.frame(switchGenes, stringsAsFactors = FALSE)
switchOut <- data.frame(
  geneId           = switchDf$geneID,
  direction        = switchDf$direction,
  fdr              = switchDf$FDR,
  pseudotime       = switchDf$switch_at_time,
  pseudoR2s        = switchDf$pseudoR2s,
  stringsAsFactors = FALSE
)

# Filter by FDR threshold
switchOut <- subset(switchOut, !is.na(fdr))
if (config$verbose) { message("Genes with valid FDR: ", nrow(switchOut)) }

switchOut <- subset(switchOut, fdr <= config$geneSwitchesFdrThreshold)
if (config$verbose) { message("Genes passing FDR threshold (", config$geneSwitchesFdrThreshold, "): ", nrow(switchOut)) }

switchOut <- switchOut[order(switchOut$fdr, switchOut$geneId, na.last = TRUE), ]

# --- Save Results ---
if (config$saveResults) {
  saveMonocle3(cds, ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
  saveObject(switchGenes, ptPaths$geneSwitches, config, "switch genes")
  saveObject(switchOut, ptPaths$geneSwitchesTsv, config, "gene switches report")
}

if (config$verbose) { message("Number of switch genes: ", nrow(switchOut)) }
message("Done!")
