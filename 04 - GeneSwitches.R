# =============================================================================
# 04 - GeneSwitches
#
# Perform logistic modeling of gene expression state transitions across
# pseudotime using the GeneSwitches package on raw data.
# Output includes filtered switch-like genes.
# =============================================================================

# --- Initialization ---
source("managers/geneSwitchesManager.R")
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
cds <- loadPseudotimeTrajectory(ptPaths, config)

# --- Prepare Raw Expression for GeneSwitches ---
if (config$verbose) { message("Preparing raw data for GeneSwitches") }
colData(cds)$Pseudotime <- pseudotime(cds)
assay(cds, "expdata")   <- log1p(assay(cds, "counts"))
expdata <- assay(cds, "expdata")

# --- Optional prefilter: detection, variance, and Spearman vs pseudotime ---
pt <- colData(cds)$Pseudotime
detect_frac <- rowMeans(expdata > 0, na.rm = TRUE)
var_ok      <- apply(expdata, 1, var, na.rm = TRUE) >= config$geneSwitchesVarThreshold

if (isTRUE(config$geneSwitchesSpearmanFilter)) {
  if (config$verbose) message("Prefiltering by |Spearman| ≥ ", config$geneSwitchesRho,
                              ", detect ≥ ", config$geneSwitchesMinDetect,
                              ", var ≥ ", config$geneSwitchesVarThreshold)
  gene_rho <- apply(expdata, 1, function(x) { suppressWarnings(cor(x, pt, method = "spearman", use = "pairwise.complete.obs")) })
  keep <- (abs(gene_rho) >= config$geneSwitchesRho) & (detect_frac >= config$geneSwitchesMinDetect) & var_ok
  if (config$verbose) message("Genes kept pre-binarization: ", sum(keep), " / ", nrow(expdata))
  if (sum(keep) == 0L) {
    warning("No genes passed the prefilter; disabling it for this run.")
  } else {
    cds     <- cds[keep, ]
    expdata <- assay(cds, "expdata")
  }
} else {
  if (config$verbose) message("Skipping Spearman prefilter; keeping genes by detection & variance only")
  keep    <- (detect_frac >= config$geneSwitchesMinDetect) & var_ok
  cds     <- cds[keep, ]
  expdata <- assay(cds, "expdata")
}

# =============================================================================
# --- Core GeneSwitches Pipeline with safety net ---
if (isTRUE(config$geneSwitchesFixedCutoff)) {

  # Fixed cutoff path
  if (config$verbose) { message("Running GeneSwitches Pipeline (fixed cutoff = ", config$geneSwitchesCutOff, ")") }
  cds <- binarize_exp(cds, fix_cutoff = TRUE, binarize_cutoff = config$geneSwitchesCutOff, ncores = config$cores)

} else {
  
  # Mixture model path (with fallback)
  if (config$verbose) message("Running GeneSwitches Pipeline (mixture model)")
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

# Fit logistic models (use correct arg name; avoid downsampling in sparse data)
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
switchOut <- subset(switchOut, !is.na(fdr))
switchOut <- switchOut[order(switchOut$fdr, switchOut$geneId, na.last = TRUE), ]

# --- Save Results ---
if (config$saveResults) {
  saveMonocle3GeneSwitches(cds, ptPaths, config)
  saveGeneSwitches(switchGenes, ptPaths, config)
  saveGeneSwitchesReport(switchOut, ptPaths, config)
}

if (config$verbose) { message("Number of switch genes: ", nrow(switchOut)) }
message("Done!")
