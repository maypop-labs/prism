# =============================================================================
# 05 - GRN.R (New Modular Version)
# Build gene regulatory network with SCENIC/GENIE3, enhanced with GeneSwitches
# =============================================================================

# --- Initialization ---
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

# --- Load Data ---
if (config$verbose) message("Loading trajectory and switch gene data...")
cds         <- loadPseudotimeTrajectory(ptPaths, config)
switchGenes <- loadSwitchGenes(ptPaths, config)

# =============================================================================
# MAIN EXECUTION PIPELINE
# =============================================================================

if (config$verbose) message("Starting enhanced GRN construction pipeline...")

# Execute modular pipeline with error handling
tryCatch({
  
  # Phase 1: Data preparation  
  prepData <- prepareGrnData(cds, switchGenes, config)
  
  # Phase 2a: GENIE3 pipeline
  scenicOptions <- runGenie3Pipeline(prepData, config)
  
  # Phase 2b: SCENIC modules and regulons  
  scenicOptions <- runScenicModules(scenicOptions, config)
  
  # Phase 2c: SCENIC scoring (needed for some metadata)
  scenicOptions <- runScenicScoring(scenicOptions, prepData$exprMatLog, config)
  
  # Phase 3: Enhanced GRN construction
  enhancedEdges <- buildEnhancedGrn(scenicOptions, prepData$switchGenes, prepData$exprMat, config)
  
  # Phase 4: Filtering and finalization
  finalGrn <- filterAndFinalizeGrn(enhancedEdges, config)
  
  # Phase 5: Output and cleanup
  saveGrnOutputs(finalGrn, ptPaths, config)
  
  if (config$verbose) message("GRN construction completed successfully!")
  
}, error = function(e) {
  message("ERROR in GRN pipeline: ", e$message)
  message("\nCheckpoint files in 'int/' directory:")
  if (dir.exists("int")) {
    intFiles <- list.files("int", pattern = "\\.Rds$")
    if (length(intFiles) > 0) {
      message("  Available: ", paste(intFiles, collapse = ", "))
      message("\nYou can resume from the last successful checkpoint by re-running the script")
    } else {
      message("  No checkpoint files found")
    }
  }
  stop(e)
})

message("Done!")
