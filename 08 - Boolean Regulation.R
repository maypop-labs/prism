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

# -- load Switch DEGs ---
if (!file.exists(ptPaths$switchDegs)) stop("Switch DEG RDS file not found: ", ptPaths$switchDegs)
message("Loading switch DEGs from: ", ptPaths$switchDegs)
switchDEGs <- readRDS(ptPaths$switchDegs)

# -- Load GRN edge list ---
if (!file.exists(ptPaths$grnPart02Edges)) stop("Edges RDS file not found: ", ptPaths$grnPart02Edges)
message("Loading edges from: ", ptPaths$grnPart02Edges)
edges <- readRDS(ptPaths$grnPart02Edges)

if (!file.exists(ptPaths$grnPart02)) stop("Final GRN RDS file not found: ", ptPaths$grnPart02)
message("Loading final GRN from: ", ptPaths$grnPart02)
g <- readRDS(ptPaths$grnPart02)

# Report
message("Total switch DEGs loaded: ", nrow(switchDEGs))
message("Total edges loaded: ", nrow(edges))

# --- Prepare binary matrix and pseudotime order ---
matBin    <- assay(cds, "binary")
cellOrder <- order(colData(cds)$Pseudotime)
nCells    <- length(cellOrder)

winSize <- max(5, floor(config$winSizePercent * nCells)) +floor(0.5 * max(5, floor(config$winSizePercent * nCells)))
kStep <- winSize

# --- Filter edges ---
if (config$grnMergeStronglyConnectedComponents) {
  gDf <- igraph::as_data_frame(g, what = "edges")
  colnames(gDf) <- c("TF", "Target", "corr", "regType")
  edges <- gDf %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = config$boolMaxRegulators) %>% ungroup()
} else {
edges <- edges %>% group_by(Target) %>% slice_max(order_by = abs(corr), n = config$boolMaxRegulators) %>% ungroup()
}

# Report
message("After edge filtering:")
message("  Edges after slice_max: ", nrow(edges))
message("  Unique targets in edges: ", length(unique(edges$Target)))
message("  Unique TFs in edges: ", length(unique(edges$TF)))


# --- Filter to genes present in expression matrix ---
geneList <- rownames(matBin)
adjTable <- edges %>% filter(TF %in% geneList & Target %in% geneList) %>% select(TF, Target)
targetGenes <- unique(adjTable$Target)

# Report
message("Genes in expression matrix: ", length(geneList))
message("After gene filtering:")
message("  Rows in adjTable: ", nrow(adjTable))
message("  Unique targets in adjTable: ", length(unique(adjTable$Target)))
message("  Unique TFs in adjTable: ", length(unique(adjTable$TF)))
message("Final target genes for Boolean inference: ", length(targetGenes))
message("  First 10 targets: ", paste(head(targetGenes, 10), collapse = ", "))

stop()

# --- Main loop for determining Boolean rules ---
message("Inferring Boolean rules for ", length(unique(adjTable$Target)), " target genes")
boolRules <- list()
geneIndex <- 1
for (gene in unique(adjTable$Target)) {
  message("[", geneIndex, "/", length(unique(adjTable$Target)), "] Gene: ", gene)
  geneIndex <- geneIndex + 1
  
  # Create regulators and immediately fix any list issues
  regulators_raw <- adjTable$TF[adjTable$Target == gene]
  
  # Force conversion to character vector
  if (is.list(regulators_raw) || class(regulators_raw)[1] == "list") {
    message("Converting regulators from list to character")
    regulators <- character(length(regulators_raw))
    for (i in seq_along(regulators_raw)) {
      if (is.list(regulators_raw[[i]])) {
        regulators[i] <- as.character(regulators_raw[[i]][[1]])
      } else {
        regulators[i] <- as.character(regulators_raw[[i]])
      }
    }
  } else {
    regulators <- as.character(regulators_raw)
  }
  regulators <- unique(regulators)

  if (length(regulators) == 0) {
    message("No regulators. Skipping.")
    next
  }
  
  # Verify regulators exist in matrix
  missing <- setdiff(regulators, rownames(matBin))
  if (length(missing) > 0) {
    message("Missing regulators: ", paste(missing, collapse = ", "), ". Skipping.")
    next
  }
  
  ioDf <- makeInputOutputPairs(
    targetGene = gene,
    regulators = regulators,
    matBin     = matBin,
    cellOrder  = cellOrder,
    k          = kStep
  )
  
  if (is.null(ioDf) || nrow(ioDf) < config$boolMinPairs) {
    message("Too few input-output pairs. Skipping.")
    next
  }
  
  # --- Find the Boolean Rules ---
  regulatorSigns <- setNames(edges$regType[edges$Target == gene], edges$TF[edges$Target == gene])
  res <- findBestBooleanRulesWithPrior(ioDf, regulatorSigns = regulatorSigns)
  pattern <- combineBooleanFunctionsByOr(res$bestFns)
  
  boolRules[[gene]] <- list(
    regulators  = regulators,
    bestScore   = res$score,
    outPattern  = pattern,
    methodUsed  = res$methodUsed,
    biologicallyPlausible = res$biologicallyPlausible
  )
}

# --- Save Rules ---
if (config$saveResults) {
  message("Saving Boolean rules to: ", ptPaths$booleanRules)
  saveRDS(boolRules, file = ptPaths$booleanRules)
  
  generateBooleanRuleReport(boolRules, edges, paths)
  ruleStats <- extractRuleStatistics(boolRules)
  
  jsonOutput <- list(
    metadata = list(
      timestamp = Sys.time(),
      parameters = config,
      summary_stats = ruleStats
    ),
    rules = boolRules,
    network_structure = edges
  )
  write_json(jsonOutput, paste0(paths$base$json, cellType, "_", trajectory, "boolean_rules_complete.json"), pretty = TRUE)
  
  # Structured TSV for tabular analysis
  rulesFlat <- data.frame(
    gene = rep(names(boolRules), sapply(boolRules, function(x) length(x$regulators))),
    regulator = unlist(sapply(boolRules, function(x) x$regulators)),
    rule_quality = rep(sapply(boolRules, function(x) x$bestScore), 
                       sapply(boolRules, function(x) length(x$regulators))),
    method_used = rep(sapply(boolRules, function(x) x$methodUsed), 
                      sapply(boolRules, function(x) length(x$regulators)))
  )
  write.table(rulesFlat, paste0(paths$base$tsv, cellType, "_", trajectory, "boolean_rules.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # NetworkX-compatible edge list
  networkxFormat <- edges %>%
    left_join(
      data.frame(
        gene = names(boolRules),
        rule_quality = sapply(boolRules, function(x) x$bestScore)
      ),
      by = c("Target" = "gene")
    )
  write.table(networkxFormat, paste0(paths$base$tsv, cellType, "_", trajectory, "network_for_ai.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
}

message("Done!")
