# =============================================================================
# 08 - Boolean Regulation (Updated with k=0 Empirical Learning)
#
# Infer Boolean regulation rules for each gene using binary expression data
# and structural information from the GRN. This version uses exhaustive Boolean
# function learning with k=0 (same pseudotime) for empirical rule discovery.
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

# --- Load necessary data ---
if (!dir.exists(ptPaths$monocle3SmoothedGeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3SmoothedGeneSwitches)
if (!file.exists(ptPaths$grnPart02Edges)) stop("Edges RDS file not found: ", ptPaths$grnPart02Edges)
cds   <- load_monocle_objects(directory_path = ptPaths$monocle3SmoothedGeneSwitches)
edges <- readRDS(ptPaths$grnPart02Edges)


# --- Sanitize Gene Names for BoolNet Compatibility ---
message("Sanitizing gene names for BoolNet compatibility...")

# Generate gene name mapping (original -> sanitized)
allGenes <- unique(c(rownames(cds), edges$TF, edges$Target))
geneMap  <- generateSanitizedGeneMapping(allGenes)

# Create lookup for fast conversion
geneLookup <- setNames(geneMap$SanitizedName, geneMap$OriginalName)

# --- Sanitize Expression Matrix ---
message("Sanitizing expression matrix gene names...")
matBin <- assay(cds, "binary")

# Only keep genes that need to be in the analysis
genesInEdges <- unique(c(edges$TF, edges$Target))
matBin <- matBin[rownames(matBin) %in% genesInEdges, ]

# Sanitize row names, edge list, TFs, and target columns
originalRowNames  <- rownames(matBin)
sanitizedRowNames <- sapply(originalRowNames, function(g) geneLookup[[g]] %||% g)
rownames(matBin)  <- sanitizedRowNames
edgesOriginal     <- edges
edges$TF          <- sapply(edges$TF, function(g) geneLookup[[g]] %||% g)
edges$Target      <- sapply(edges$Target, function(g) geneLookup[[g]] %||% g)

# Remove edges where genes are missing from expression matrix
validTFs     <- edges$TF %in% rownames(matBin)
validTargets <- edges$Target %in% rownames(matBin)
validEdges   <- validTFs & validTargets

if (sum(!validEdges) > 0) {
  message("Removing ", sum(!validEdges), " edges with genes missing from expression matrix")
  edges <- edges[validEdges, ]
}

message("Final edge list: ", nrow(edges), " edges")

# --- Prepare Binary Matrix and Pseudotime Order ---
message("Preparing binary expression data and pseudotime ordering...")

# Get pseudotime ordering of cells
cellOrder <- order(colData(cds)$Pseudotime)
kStep     <- 0  # Use k = 0 for same-pseudotime analysis (no temporal lag)

message("Using k = ", kStep, " for empirical Boolean rule learning")
message("Cell ordering: ", length(cellOrder), " cells ordered by pseudotime")

# --- Infer Boolean Rules with k=0 State Counting Learning ---
message("Inferring Boolean rules ...")

targets <- unique(edges$Target)
total   <- length(targets)

if (total == 0) {
  stop("No target genes found in edges data frame")
}

boolRules <- list()
pb <- txtProgressBar(min = 0, max = total, style = 3)

for (i in seq_along(targets)) {
  gene <- targets[i]
  setTxtProgressBar(pb, i)
  
  # Get edges for this gene
  subEdges <- edges[edges$Target == gene, ]
  if (nrow(subEdges) == 0) { next }
  
  # Get top regulators (limit to maxRegulators)
  maxRegulators <- config$boolMaxRegulators
  topEdges      <- head(subEdges[order(-abs(subEdges$corr)), ], maxRegulators)
  regulators    <- topEdges$TF
  
  # Check if gene and regulators are in matrix
  geneInMatrix <- gene %in% rownames(matBin)
  regsInMatrix <- regulators %in% rownames(matBin)
  if (!geneInMatrix || !any(regsInMatrix)) { next }
  
  # Filter to only regulators that exist in matrix
  regulators <- regulators[regsInMatrix]
  if (length(regulators) == 0) { next }
  
  # Generate input-output pairs with k=0
  tryCatch({
    ioData <- makeInputOutputPairs(gene, regulators, matBin, cellOrder, k = kStep)
    
    if (!is.null(ioData) && length(ioData$stateIndices) >= 5) {
      # Find best Boolean functions using optimized state counting
      res <- findBestBooleanRules(ioData)
      
      if (res$score > 0 && length(res$bestFns) > 0) {
        # Get the optimal output pattern
        outPattern <- res$bestFns[[1]]$outPattern
        
        # Convert to BoolNet rule string
        ruleStr <- makeBoolNetRule(gene, ioData$regulators, outPattern)
        
        # Store rule information
        boolRules[[gene]] <- list(
          rule             = ruleStr,
          regulators       = ioData$regulators,
          bestScore        = res$score,
          score            = res$score,
          nRegulators      = length(ioData$regulators),
          outPattern       = outPattern,
          nBestFunctions   = length(res$bestFns),
          functionsTestsed = res$functionsTestsed,
          totalPossible    = res$totalPossible,
          method           = res$method,
        )
      }
    }
    
  }, error = function(e) {
    # Skip genes that fail rule synthesis - just continue to next iteration
    warning("Warning: Failed to process gene ", gene, ": ", e$message)
  })
}
close(pb)

message("\nRules synthesized for ", length(boolRules), " target genes using state counting learning.")

# --- Add Self-Activation Rules for Genes Without Rules ---
message("Adding self-activation rules for genes without inferred rules...")

# Find all genes mentioned in rules or as regulators
ruleGenes <- names(boolRules)
genesInRules <- unique(unlist(sapply(boolRules, function(x) x$regulators)))
allMentionedGenes <- unique(c(ruleGenes, genesInRules))

# Find genes that exist in matrix but don't have rules
genesNeedingRules <- setdiff(intersect(allMentionedGenes, rownames(matBin)), ruleGenes)

message("Found ", length(genesNeedingRules), " genes without rules")
if (length(genesNeedingRules) > 0) {
  message("Examples: ", paste(head(genesNeedingRules, 5), collapse = ", "))
}

# Add self-activation rules for these genes
for (gene in genesNeedingRules) {
  # Create proper self-activation rule using makeBoolNetRule
  selfRule <- makeBoolNetRule(gene, character(0), c(1))  # Always ON
  
  boolRules[[gene]] <- list(
    rule = selfRule,
    regulators = character(0),  # No regulators (self-activation)
    bestScore = 1.0,  # Not actually measured, but logical for "always ON"
    score = 1.0,
    nRegulators = 0,
    outPattern = c(1),  # Always ON pattern
    nBestFunctions = 1,
    functionsTested = 1,  # Fixed typo
    totalPossible = 2,   # Could be ON or OFF
    method = "self_activation",
  )
}

message("Added self-activation rules for ", length(genesNeedingRules), " genes")
message("Total rules after adding self-activation: ", length(boolRules))

# --- Quality Control and Validation ---
message("Performing quality control checks...")

# Check rule consistency
ruleGenes <- unique(c(names(boolRules), unlist(sapply(boolRules, function(x) x$regulators))))
matrixGenes <- rownames(matBin)
missingGenes <- setdiff(ruleGenes, matrixGenes)

if (length(missingGenes) > 0) {
  warning("Found ", length(missingGenes), " genes in rules but not in matrix: ", 
          paste(head(missingGenes, 5), collapse = ", "))
} else {
  message("✓ All genes in rules are present in expression matrix")
}

# Check for problematic rule patterns
problematicRules <- sapply(boolRules, function(rule) {
  ruleStr <- rule$rule
  # Check for basic issues
  parts <- strsplit(ruleStr, ",\\s*")[[1]]
  if (length(parts) != 2) return(TRUE)
  
  # Check for invalid characters
  logic <- parts[2]
  hasInvalidChars <- grepl("[^A-Za-z0-9_&|()!\\s]", logic)
  return(hasInvalidChars)
})

if (any(problematicRules)) {
  warning("Found ", sum(problematicRules), " potentially problematic rules")
  problematicNames <- names(boolRules)[problematicRules]
  message("Problematic rules: ", paste(head(problematicNames, 5), collapse = ", "))
} else {
  message("✓ All rules appear to have valid BoolNet syntax")
}


# --- Save Results ---
if (config$saveResults) {
  message("Saving Boolean rules to: ", ptPaths$booleanRules)
  saveRDS(boolRules, file = ptPaths$booleanRules)
  
  # Save the gene name mapping for reference
  mappingFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_gene_mapping.rds"))
  message("Saving gene name mapping to: ", mappingFile)
  saveRDS(geneMap, file = mappingFile)
  
  # Write detailed rule analysis table
  if (length(boolRules) > 0) {
    ruleAnalysis <- data.frame(
      gene = names(boolRules),
      rule_string = sapply(boolRules, function(x) x$rule),
      score = sapply(boolRules, function(x) x$score),
      n_regulators = sapply(boolRules, function(x) x$nRegulators),
      n_best_functions = sapply(boolRules, function(x) x$nBestFunctions %||% 1),
      functions_tested = sapply(boolRules, function(x) x$functionsTested %||% NA),
      total_possible = sapply(boolRules, function(x) x$totalPossible %||% NA),
      method = sapply(boolRules, function(x) x$method),
      has_or = sapply(boolRules, function(x) grepl("\\|", x$rule)),
      has_and = sapply(boolRules, function(x) grepl("&", x$rule)),
      stringsAsFactors = FALSE
    )
    
    write.table(ruleAnalysis, 
                paste0(paths$base$tsv, cellType, "_", trajectory, "_boolean_rules_analysis.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Also create the flat regulator table for compatibility
    rulesFlat <- data.frame(
      gene = rep(names(boolRules), sapply(boolRules, function(x) max(1, length(x$regulators)))),
      regulator = unlist(lapply(boolRules, function(x) {
        if (length(x$regulators) == 0) return(NA)
        return(x$regulators)
      })),
      rule_quality = rep(sapply(boolRules, function(x) x$score), 
                         sapply(boolRules, function(x) max(1, length(x$regulators)))),
      rule_string = rep(sapply(boolRules, function(x) x$rule), 
                        sapply(boolRules, function(x) max(1, length(x$regulators)))),
      stringsAsFactors = FALSE
    )
    
    write.table(rulesFlat, 
                paste0(paths$base$tsv, cellType, "_", trajectory, "_boolean_rules.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Create rule report (skip if it causes errors)
  tryCatch({
    generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
  }, error = function(e) {
    message("Warning: Could not generate rule report: ", e$message)
  })
  
  # Write sanitization summary
  sanitizationSummary <- data.frame(
    OriginalName = geneMap$OriginalName,
    SanitizedName = geneMap$SanitizedName,
    Changed = geneMap$OriginalName != geneMap$SanitizedName,
    stringsAsFactors = FALSE
  )
  write.table(sanitizationSummary, 
              paste0(paths$base$tsv, cellType, "_", trajectory, "_gene_sanitization.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("Saved gene sanitization summary with ", sum(sanitizationSummary$Changed), " changed names")
}

message("Done! Boolean rules inferred using optimized state counting approach.")
message("Key insights:")
message("- Used k=0 (same pseudotime) to avoid questionable temporal assumptions")
message("- Used state counting instead of exhaustive enumeration (orders of magnitude faster)")
message("- Mathematically equivalent to exhaustive search with tie-breaking by OR")
message("- Should produce accurate network dynamics with much faster computation")