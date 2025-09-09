# =============================================================================
# 07 - BoolNet (Boolean Network Analysis)
# Purpose: Perform Boolean network analysis and attractor identification
# =============================================================================

# Core Philosophy:
# - Load validated Boolean rules from script 06
# - Run BoolNet attractor analysis on the entire regulatory network
# - Use configurable thresholds for exhaustive vs sampling approaches
# - Preserve all quality control and validation from upstream processing

# --- Initialization ---
# source("managers/booleanManager.R")
# source("managers/booleanReportManager.R")
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

# === STAGE 1: Load Validated Data ===
message("=== STAGE 1: Loading validated Boolean rules and trajectory data ===")
cds       <- loadPseudotimeTrajectory(ptPaths, config)
boolRules <- loadBooleanRules(ptPaths, config)

message("Data loaded successfully:")
message("  - CDS object: ", nrow(cds), " genes, ", ncol(cds), " cells")
message("  - Boolean rules: ", length(boolRules), " validated rules from script 06")

# === STAGE 2: Boolean Network Analysis ===
message("\n=== STAGE 2: Creating BoolNet network and finding attractors ===")

# Generate BoolNet-compatible rule table from validated rules
message("Generating BoolNet rule table...")
ruleLines <- c("targets,factors")
for (gene in names(boolRules)) {
  ruleLines <- c(ruleLines, boolRules[[gene]]$rule)
}

message("BoolNet rule table: ", length(ruleLines) - 1, " rules for ", length(boolRules), " genes")

# Create BoolNet network object
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)
boolNetwork <- loadNetwork(tempFile)
boolNetwork$type <- "synchronous"

nNetworkGenes <- length(boolNetwork$genes)
message("BoolNet network created: ", nNetworkGenes, " genes")

# Choose analysis method based on network size and configuration
if (nNetworkGenes <= config$boolExhaustiveLimit) {
  message("Using exhaustive search (≤", config$boolExhaustiveLimit, " genes)")
  attractors <- getAttractors(
    boolNetwork,
    method      = "exhaustive",
    type        = "synchronous",
    returnTable = TRUE
  )
} else { 
  adaptiveSamples <- min(config$boolMaxSamples, max(50, nNetworkGenes * 4))
  
  message("Using random sampling: ", adaptiveSamples, " samples")
  attractors <- getAttractors(
    boolNetwork,
    method      = "random",
    startStates = adaptiveSamples,
    type        = "synchronous",
    returnTable = TRUE
  )
}

# === STAGE 3: Analyze and Report Results ===
message("\n=== STAGE 3: Analyzing attractor results ===")

nAttractors <- length(attractors$attractors)
attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))

message("✓ SUCCESS: Identified ", nAttractors, " attractors")
message("  Attractor sizes: ", paste(attractorSizes, collapse = ", "))

# Additional analysis
if (nAttractors > 0) {
  maxAttractorSize <- max(attractorSizes)
  avgAttractorSize <- round(mean(attractorSizes), 2)
  
  message("  Maximum attractor size: ", maxAttractorSize)
  message("  Average attractor size: ", avgAttractorSize)
  
  # Point attractors (steady states)
  pointAttractors <- sum(attractorSizes == 1)
  if (pointAttractors > 0) {
    message("  Point attractors (steady states): ", pointAttractors)
  }
  
  # Cycle attractors
  cycleAttractors <- sum(attractorSizes > 1)
  if (cycleAttractors > 0) {
    message("  Cycle attractors: ", cycleAttractors)
  }
}

# === STAGE 4: Save Results ===
if (config$saveResults) {
  message("\n=== STAGE 4: Saving results ===")
  
  saveAttractors(attractors, ptPaths, config)
  saveBoolNetwork(boolNetwork, ptPaths, config)
  
  # Generate analysis summary
  boolNetReport <- data.frame(
    TotalGenes = nNetworkGenes,
    TotalAttractors = nAttractors,
    PointAttractors = sum(attractorSizes == 1),
    CycleAttractors = sum(attractorSizes > 1),
    MaxAttractorSize = ifelse(nAttractors > 0, max(attractorSizes), 0),
    AvgAttractorSize = ifelse(nAttractors > 0, round(mean(attractorSizes), 2), 0),
    AnalysisMethod = ifelse(nNetworkGenes <= config$boolExhaustiveLimit, "Exhaustive", "Sampling"),
    stringsAsFactors = FALSE
  )
  saveBoolNetReport(boolNetReport, ptPaths, config)
}

# === FINAL SUMMARY ===
message("\n", paste(rep("=", 60), collapse = ""))
message("BOOLEAN NETWORK ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("SUCCESS: Analyzed Boolean network with ", nNetworkGenes, " genes")
message("Identified ", nAttractors, " attractors using ", 
        ifelse(nNetworkGenes <= config$boolExhaustiveLimit, "exhaustive", "sampling"), " approach")
message("\n", paste(rep("=", 60), collapse = ""))
message("Done!")
