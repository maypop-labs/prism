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

# --- Generate Sanitized Gene Mapping ---
geneMap <- generateSanitizedGeneMapping(rownames(cds))

# --- Filter Rules by Score ---
message("Starting with ", length(boolRules), " Boolean rules")
minScore <- config$boolNetMinScore %||% 0.65
message("Filtering rules with score ≥ ", minScore)
boolRules <- Filter(function(rule) {
  !is.null(rule$score) && rule$score >= minScore
}, boolRules)
message("Retained ", length(boolRules), " high-confidence Boolean rules")



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

# Progress bar for rule table building
geneNames <- names(boolRules)
total <- length(geneNames)
pb <- txtProgressBar(min = 0, max = total, style = 3)

for (i in seq_along(geneNames)) {
  gene <- geneNames[i]
  setTxtProgressBar(pb, i)
  
  ruleStr <- sanitizeRule(boolRules[[gene]]$rule, geneMap)
  if (grepl(",,", ruleStr) || grepl("\\bNA\\b", ruleStr)) {
    warning("Skipping malformed rule for gene: ", gene)
    next
  }
  ruleLines <- c(ruleLines, ruleStr)
}

close(pb)
message("BoolNet rule table complete: ", length(ruleLines) - 1, " rules")

# --- Convert to BoolNet Object ---
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)

tryCatch({
  boolnet <- loadNetwork(tempFile)
  message("BoolNet object created successfully")
  boolnet$type <- "synchronous"
}, error = function(e) {
  message("FAILED to create BoolNet object: ", e$message)
  stop("BoolNet creation failed: ", e$message)
})

# --- Compute Attractors ---



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