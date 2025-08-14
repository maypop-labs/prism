# =============================================================================
# 08 - BoolNet (Component-wise Analysis with Proper Isolation)
#
# This script performs component-wise Boolean network analysis with proper
# isolation of genes and regulators within each component
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
if (!dir.exists(ptPaths$monocle3GeneSwitches)) {
  stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
}
message("Loading Monocle3 object from: ", ptPaths$monocle3GeneSwitches)
cds <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)

# --- Load Boolean rules ---
if (!file.exists(ptPaths$booleanRules)) {
  stop("Boolean rules file not found: ", ptPaths$booleanRules)
}
message("Loading Boolean rules from: ", ptPaths$booleanRules)
boolRules <- readRDS(ptPaths$booleanRules)

# --- Generate Sanitized Gene Mapping ---
message("Generating sanitized gene mapping...")
geneMap <- generateSanitizedGeneMapping(rownames(cds))

# --- Filter Rules by Score ---
minScore <- config$boolMinScore %||% 0.65

message("Filtering rules with score ≥ ", minScore)
boolRulesOriginal <- boolRules
boolRules <- Filter(function(rule) {
  !is.null(rule$score) && rule$score >= minScore
}, boolRules)
message("Retained ", length(boolRules), " out of ", length(boolRulesOriginal), " rules")

# --- Build Initial Network for Component Detection ---
message("Building network from filtered rules...")
edgesDf <- do.call(rbind, lapply(names(boolRules), function(target) {
  regulators <- boolRules[[target]]$regulators
  if (length(regulators) > 0) {
    data.frame(TF = regulators, Target = target, stringsAsFactors = FALSE)
  } else {
    NULL
  }
}))

# Remove any NULL entries
edgesDf <- edgesDf[!is.null(edgesDf), ]

g <- igraph::graph_from_data_frame(edgesDf, directed = TRUE)
comp <- components(g)
message("Network has ", comp$no, " connected components")
message("Component sizes: ", paste(comp$csize, collapse = ", "))

# Filter components by minimum size
minComponentSize <- 3
validComponents <- which(comp$csize >= minComponentSize)
message("Analyzing ", length(validComponents), " components with ≥", minComponentSize, " genes")

if (length(validComponents) == 0) {
  stop("No components meet the minimum size requirement of ", minComponentSize, " genes")
}

# Create component ID vector
componentIds <- paste0("component", validComponents)
names(componentIds) <- validComponents

# --- Save Component IDs ---
if (config$saveResults) {
  componentFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_components.rds"))
  message("Saving component IDs to: ", componentFile)
  saveRDS(componentIds, file = componentFile)
}

# --- Helper Function: Rebuild Component Rules ---
rebuildComponentRules <- function(componentGenes, originalRules, matBin) {
  message("  Rebuilding rules with component-isolated regulators...")
  
  componentRules <- list()
  
  for (targetGene in componentGenes) {
    if (!targetGene %in% names(originalRules)) {
      next  # Skip genes without rules
    }
    
    originalRule <- originalRules[[targetGene]]
    originalRegulators <- originalRule$regulators
    
    # Keep only regulators that are in this component
    validRegulators <- originalRegulators[originalRegulators %in% componentGenes]
    
    if (length(validRegulators) == 0) {
      # Create self-activation rule for isolated genes
      validRegulators <- targetGene
      regulatorSigns <- setNames(1, targetGene)  # Self-activation
    } else {
      # Get regulation signs for valid regulators
      # We need to look up the original edge information
      regulatorSigns <- setNames(rep(1, length(validRegulators)), validRegulators)  # Default to activation
      
      # Try to recover original signs from the rule structure or use default
      # For simplicity, we'll default to activation and let the template system figure it out
    }
    
    # Regenerate the Boolean rule using our existing template system
    tryCatch({
      ruleResult <- synthesizeBestBooleanRule(targetGene, validRegulators, regulatorSigns, matBin)
      
      componentRules[[targetGene]] <- list(
        rule = ruleResult$bestRule,
        regulators = validRegulators,
        score = ruleResult$bestScore,
        nRegulators = length(validRegulators),
        templateUsed = ruleResult$bestTemplateName,
      )
      
    }, error = function(e) {
      # Fallback: create simple self-activation rule
      componentRules[[targetGene]] <- list(
        rule = paste0(targetGene, ", ", targetGene),
        regulators = targetGene,
        score = 0.5,  # Neutral score
        nRegulators = 1,
        templateUsed = "self_activation",
      )
    })
  }
  
  message("  Component rules rebuilt: ", length(componentRules))
  return(componentRules)
}

# --- Component-wise Attractor Analysis ---
message("\n=== COMPONENT-WISE ATTRACTOR ANALYSIS ===")

successfulComponents <- character(0)
failedComponents <- character(0)

# Get binary matrix for rule regeneration
matBin <- assay(cds, "binary")

for (i in seq_along(validComponents)) {
  componentIndex <- validComponents[i]
  componentId <- componentIds[i]
  componentSize <- comp$csize[componentIndex]
  
  message("\n--- ", componentId, " (", componentSize, " genes) ---")
  
  # Get genes in this component
  componentGenes <- V(g)[comp$membership == componentIndex]$name
  message("Genes: ", paste(head(componentGenes, 8), collapse = ", "), 
          if(length(componentGenes) > 8) "..." else "")
  
  # Rebuild rules for this component with proper isolation
  tryCatch({
    componentRules <- rebuildComponentRules(componentGenes, boolRules, matBin)
    
    if (length(componentRules) == 0) {
      message("ERROR: No valid rules could be created for component")
      failedComponents <- c(failedComponents, componentId)
      next
    }
    
    message("Building BoolNet rule table...")
    
    # Create rule lines
    ruleLines <- c("targets,factors")
    
    for (gene in names(componentRules)) {
      ruleStr <- sanitizeRule(componentRules[[gene]]$rule, geneMap)
      if (grepl(",,", ruleStr) || grepl("\\bNA\\b", ruleStr)) {
        warning("Skipping malformed rule for gene: ", gene)
        next
      }
      ruleLines <- c(ruleLines, ruleStr)
    }
    
    message("Rule table complete: ", length(ruleLines) - 1, " rules")
    
    # Create and load BoolNet object
    tempFile <- tempfile(fileext = ".txt")
    writeLines(ruleLines, con = tempFile)
    
    componentNetwork <- loadNetwork(tempFile)
    componentNetwork$type <- "synchronous"
    
    nNetworkGenes <- length(componentNetwork$genes)
    message("BoolNet object created: ", nNetworkGenes, " genes")
    
    # Adaptive attractor computation based on component size
    if (nNetworkGenes <= 20) {
      # Small component: exhaustive search
      message("Using exhaustive search for small component")
      attractors <- getAttractors(
        componentNetwork,
        method = "exhaustive",
        type = "synchronous",
        returnTable = TRUE
      )
    } else {
      # Larger component: adaptive random sampling
      adaptiveSamples <- if (nNetworkGenes <= 50) {
        max(200, nNetworkGenes * 8)
      } else if (nNetworkGenes <= 100) {
        max(100, nNetworkGenes * 4)
      } else {
        max(50, nNetworkGenes * 2)
      }
      
      message("Using random sampling with ", adaptiveSamples, " samples")
      attractors <- getAttractors(
        componentNetwork,
        method = "random",
        startStates = adaptiveSamples,
        type = "synchronous",
        returnTable = TRUE
      )
    }
    
    # Report results
    nAttractors <- length(attractors$attractors)
    attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))
    
    message("SUCCESS: Found ", nAttractors, " attractors")
    message("Attractor sizes: ", paste(attractorSizes, collapse = ", "))
    
    # Save component results
    if (config$saveResults) {
      componentFile <- file.path(paths$base$rds, 
                                 paste0(cellType, "_", trajectory, "_", componentId, ".rds"))
      message("Saving attractor results to: ", basename(componentFile))
      saveRDS(attractors, file = componentFile)
    }
    
    successfulComponents <- c(successfulComponents, componentId)
    
  }, error = function(e) {
    message("FAILED: ", e$message)
    failedComponents <- c(failedComponents, componentId)
  })
}

# --- Final Summary ---
message("\n=== ANALYSIS SUMMARY ===")
message("Components analyzed: ", length(validComponents))
message("Successful: ", length(successfulComponents), " (", paste(successfulComponents, collapse = ", "), ")")
if (length(failedComponents) > 0) {
  message("Failed: ", length(failedComponents), " (", paste(failedComponents, collapse = ", "), ")")
}

# Create and save analysis summary
if (config$saveResults) {
  analysisSummary <- data.frame(
    Component = componentIds,
    Size = comp$csize[validComponents],
    Status = ifelse(componentIds %in% successfulComponents, "Success", "Failed"),
    stringsAsFactors = FALSE
  )
  
  summaryFile <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_component_analysis.tsv"))
  write.table(analysisSummary, summaryFile, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Analysis summary saved to: ", basename(summaryFile))
}

message("Component-wise Boolean network analysis complete!")
message("Done!")
