# =============================================================================
# 09 - BoolNet (Refactored with Connectivity Check and Robust Sanitization)
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

# --- Load Boolean rules ---
if (!file.exists(ptPaths$booleanRules)) stop("Boolean rules file not found: ", ptPaths$booleanRules)
message("Loading Boolean rules from: ", ptPaths$booleanRules)
boolRules <- readRDS(ptPaths$booleanRules)

# --- Filter Rules by Score ---
minScore <- config$boolNetMinScore %||% 0.65
message("Filtering rules with score ≥ ", minScore)
boolRules <- Filter(function(rule) {
  !is.null(rule$score) && rule$score >= minScore
}, boolRules)
message("Retained ", length(boolRules), " high-confidence rules")

# --- Generate Sanitized Gene Mapping ---
geneMap <- generateSanitizedGeneMapping(rownames(cds))

# --- Check Connectivity of Filtered Rules ---
edgesDf <- do.call(rbind, lapply(names(boolRules), function(target) {
  data.frame(TF = boolRules[[target]]$regulators, Target = target, stringsAsFactors = FALSE)
}))
g <- igraph::graph_from_data_frame(edgesDf, directed = TRUE)
comp <- components(g)
message("Network has ", comp$no, " connected components")

if (comp$no > 1) {
  largest <- which.max(comp$csize)
  disconnected <- V(g)[comp$membership != largest]
  message("Warning: ", length(disconnected), " nodes are outside the largest component")
}

# --- Build BoolNet Rule Table ---
message("Building BoolNet rule table")
ruleLines <- c("targets,factors")
for (gene in names(boolRules)) {
  ruleStr <- sanitizeRule(boolRules[[gene]]$rule, geneMap)
  if (grepl(",,", ruleStr) || grepl("\\bNA\\b", ruleStr)) {
    warning("Skipping malformed rule for gene: ", gene)
    next
  }
  ruleLines <- c(ruleLines, ruleStr)
}
message("BoolNet rule table complete: ", length(ruleLines) - 1, " rules")

# ---

# --- DEBUG: Inspect rule table ---
message("=== DEBUGGING RULE TABLE ===")
message("First 10 rule lines:")
for (i in 1:min(10, length(ruleLines))) {
  message(i, ": ", ruleLines[i])
}

# Check for common issues
# Look for actual problems, not just gene names containing "NA"
problematicRules <- ruleLines[grepl("\\bNA\\b|,,|^\\s*$|^targets,factors$", ruleLines)]
problematicRules <- problematicRules[!grepl("^targets,factors$", problematicRules)]  # Exclude header
if (length(problematicRules) > 0) {
  message("Found problematic rules:")
  for (rule in problematicRules) {
    message("  PROBLEM: ", rule)
  }
}

# Check BoolNet object
message("BoolNet genes: ", length(boolnet$genes))
message("BoolNet genes: ", paste(head(boolnet$genes, 10), collapse = ", "))
message("==============================")

# ---

# --- Convert to BoolNet Object ---
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)

# DEBUG: Show the temp file content
message("=== TEMP FILE CONTENT ===")
message("File location: ", tempFile)
message("First 10 lines of temp file:")
tempContent <- readLines(tempFile)
for (i in 1:min(10, length(tempContent))) {
  message(i, ": ", tempContent[i])
}
message("==========================")

# Try to load with error handling
tryCatch({
  boolnet <- loadNetwork(tempFile)
  message("BoolNet object created successfully")
  boolnet$type <- "synchronous"
}, error = function(e) {
  message("FAILED to create BoolNet object: ", e$message)
  message("Temp file still available at: ", tempFile)
  stop("BoolNet creation failed")
})

# --- Compute Attractors ---
message("Computing attractors with BoolNet")
tryCatch({
  attractors <- getAttractors(
    boolnet,
    method = "random",
    startStates = config$boolNetStartStates,
    returnTable = TRUE,
    type = "synchronous"
  )
}, error = function(e) {
  stop("Attractor computation failed: ", e$message)
})

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