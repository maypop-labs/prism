# =============================================================================
# 05 - GRN.R
# Build gene regulatory network with SCENIC/GENIE3, integrated with GeneSwitches
#
# Note: Possible > 14-hour runtime for new pseudotime trajectories.
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

# --- SCENIC Recovery Logic ---
if (!config$grnFromScratch) {
  # Check if int and output exist locally
  localIntExists    <- dir.exists("int")
  localOutputExists <- dir.exists("output")
  
  if (!localIntExists || !localOutputExists) {
    # Check if they exist in trajectory-specific scenic folder
    scenicIntPath    <- file.path(ptPaths$scenic, "int")
    scenicOutputPath <- file.path(ptPaths$scenic, "output")
    
    if (dir.exists(scenicIntPath) || dir.exists(scenicOutputPath)) {
      # Copy from scenic folder to local
      if (dir.exists(scenicIntPath) && !localIntExists) {
        file.copy(scenicIntPath, ".", recursive = TRUE)
      }
      if (dir.exists(scenicOutputPath) && !localOutputExists) {
        file.copy(scenicOutputPath, ".", recursive = TRUE)
      }
      if (config$verbose) message("Restored SCENIC files from previous run")
    }
  }
}

# --- Load Data ---
cds         <- loadMonocle3(ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory with binary assay")
switchGenes <- loadObject(ptPaths$geneSwitches, config, "switch genes")

# =============================================================================
# MAIN EXECUTION PIPELINE
# =============================================================================

if (config$verbose) message("Starting GRN construction pipeline...")

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
  
  # Phase 3: GRN construction (now uses binary correlations)
  # Extract binary matrix from CDS for correlation analysis
  matBin <- assay(cds, "binary")
  grnEdges <- buildGrn(scenicOptions, prepData$switchGenes, prepData$exprMat, matBin, config)
  
  # Phase 4: Filtering and finalization
  finalGrn <- filterAndFinalizeGrn(grnEdges, config)
  
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
