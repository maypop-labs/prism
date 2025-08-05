# =============================================================================
# 08 - Boolean Regulation
#
# Infer Boolean regulation rules for each gene using binary expression data
# and structural information from the GRN. Output a list of Boolean rules.
#
# This script has a long runtime. Grab a cup of coffee. C[_]
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

# -- Load GRN edge list ---
if (!file.exists(ptPaths$grnPart02Edges)) stop("Edges RDS file not found: ", ptPaths$grnPart02Edges)
message("Loading edges from: ", ptPaths$grnPart02Edges)
edges <- readRDS(ptPaths$grnPart02Edges)

# --- Synthesize Boolean Rules ---
message("Synthesizing Boolean rules from SCENIC edges...")
matBin    <- assay(cds, "binary")
boolRules <- synthesizeBooleanRulesBatch(
  edges = edges,
  matBin = matBin,
  maxRegulators = config$boolMaxRegulators %||% 5
)

message("Rules synthesized for ", length(boolRules), " target genes.")

# --- Save Rules ---
if (config$saveResults) {
  message("Saving Boolean rules to: ", ptPaths$booleanRules)
  saveRDS(boolRules, file = ptPaths$booleanRules)
  
  # Write flat TSV table
  if (length(boolRules) > 0) {
    rulesFlat <- data.frame(
      gene = rep(names(boolRules), sapply(boolRules, function(x) length(x$regulators))),
      regulator = unlist(sapply(boolRules, function(x) x$regulators)),
      rule_quality = rep(sapply(boolRules, function(x) x$score), sapply(boolRules, function(x) length(x$regulators)))
    )
    write.table(rulesFlat, paste0(paths$base$tsv, cellType, "_", trajectory, "_boolean_rules.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Create rule report
  generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
}

message("Done!")
