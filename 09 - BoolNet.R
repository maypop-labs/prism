# =============================================================================
# 09 - BoolNet
#
# Convert Boolean rules into a BoolNet object, compute attractors,
# and save the BoolNet network, attractors, and sanitized gene name mapping.
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/booleanReportManager.R")
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

if (!file.exists(ptPaths$booleanRules)) stop("Boolean rules file not found: ", ptPaths$booleanRules)
message("Loading Monocle3 object from: ", ptPaths$monocle3SmoothedGeneSwitches)
boolRules <- readRDS(ptPaths$booleanRules)

geneMap <- generateSanitizedGeneMapping(rownames(cds))

# --- Build BoolNet Rule Table ---
message("Building BoolNet rule table")
ruleLines <- c("targets,factors")
for (gene in names(boolRules)) {
  ruleStr <- sanitizeRule(boolRules[[gene]]$rule, geneMap)
  ruleLines <- c(ruleLines, ruleStr)
}

# --- Convert to BoolNet Object ---
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)
boolnet <- loadNetwork(tempFile)
boolnet$type <- "synchronous"

# --- Compute Attractors ---
message("Computing attractors with BoolNet")
attractors <- getAttractors(
  boolnet,
  method = "random",
  startStates = config$boolNetStartStates,
  returnTable = TRUE,
  type = "synchronous"
)

# --- Save Outputs ---
if (config$saveResults) {
  message("Saving BoolNet object to: ", ptPaths$boolnet)
  saveRDS(boolnet, file = ptPaths$boolnet)

  message("Saving attractors to: ", ptPaths$attractors)
  saveRDS(attractors, file = ptPaths$attractors)

  message("Saving gene mapping to: ", ptPaths$geneMap)
  saveRDS(geneMap, file = ptPaths$geneMap)
}

message("Done!")
