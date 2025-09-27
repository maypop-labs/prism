# =============================================================================
# 06 - Boolean Regulation
# Purpose: Generate production-quality Boolean rules using `i + k` + SCENIC data
# =============================================================================

source("managers/booleanManager.R")
source("managers/booleanReportManager.R")
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

# --- STAGE 1: Robust Data Loading and Validation ---
cds   <- loadMonocle3(ptPaths$monocle3GeneSwitches, config, "GeneSwitches trajectory")
edges <- loadObject(ptPaths$grnEdges, config, "GRN edges")

# Validate SCENIC metadata columns
required_cols <- c("TF", "Target", "corr")
missing_cols <- setdiff(required_cols, colnames(edges))
if (length(missing_cols) > 0) {
  stop("Missing required columns in edges: ", paste(missing_cols, collapse = ", "))
}

# Add missing optional columns with defaults
if (!"regType" %in% colnames(edges)) {
  message("Adding default regType (Activation) - consider running SCENIC with regulatory direction")
  edges$regType <- "Activation"
}
if (!"motifConfidence" %in% colnames(edges)) { edges$motifConfidence <- NA }
if (!"NES" %in% colnames(edges)) { edges$NES <- NA }
if (!"hasMotif" %in% colnames(edges)) { edges$hasMotif <- NA }

message("Data validation complete:")
message("  - CDS object: ", nrow(cds), " genes, ", ncol(cds), " cells")
message("  - Edge list: ", nrow(edges), " regulatory relationships")
message("  - SCENIC columns: ", paste(colnames(edges), collapse = ", "))

# --- STAGE 2: Comprehensive Gene Name Sanitization ---
message("=== STAGE 2: Sanitizing gene names for BoolNet compatibility ===")

# Generate comprehensive gene name mapping
allGenes   <- unique(c(rownames(cds), edges$TF, edges$Target))
geneMap    <- generateSanitizedGeneMapping(allGenes)
geneLookup <- setNames(geneMap$SanitizedName, geneMap$OriginalName)

message("Gene name mapping:")
message("  - Total unique genes: ", length(allGenes))
message("  - Genes requiring sanitization: ", sum(geneMap$OriginalName != geneMap$SanitizedName))

# Apply sanitization to expression matrix
matBin            <- assay(cds, "binary")
genesInEdges      <- unique(c(edges$TF, edges$Target))
matBin            <- matBin[rownames(matBin) %in% genesInEdges, ]
originalRowNames  <- rownames(matBin)
sanitizedRowNames <- sapply(originalRowNames, function(g) geneLookup[[g]] %||% g)
rownames(matBin)  <- sanitizedRowNames
edges$TF          <- sapply(edges$TF, function(g) geneLookup[[g]] %||% g)
edges$Target      <- sapply(edges$Target, function(g) geneLookup[[g]] %||% g)
validTFs          <- edges$TF %in% rownames(matBin)
validTargets      <- edges$Target %in% rownames(matBin)
validEdges        <- validTFs & validTargets

if (sum(!validEdges) > 0) {
  message("Removing ", sum(!validEdges), " edges with genes missing from expression matrix")
  edges <- edges[validEdges, ]
}

message("Final dataset:")
message("  - Expression matrix: ", nrow(matBin), " genes, ", ncol(matBin), " cells")
message("  - Valid edges: ", nrow(edges), " regulatory relationships")

# --- STAGE 3: Prepare Binary Data and Cell Ordering ---
message("=== STAGE 3: Preparing binary expression data and pseudotime ordering ===")

cellOrder <- order(colData(cds)$Pseudotime)

message("Configuration:")
message("  - Using k = ", config$boolKValue)
message("  - Cell ordering: ", length(cellOrder), " cells ordered by pseudotime")
message("  - Max regulators per gene: ", config$boolMaxRegulators)
message("  - Minimum rule score: ", config$boolMinScore)

# --- STAGE 4: Boolean Rule Inference ---
message("=== STAGE 4: SCENIC-integrated Boolean rule inference ===")

# FIXED: Create rules for all genes mentioned in network (TFs + Targets) to ensure complete coverage
allNodes <- union(unique(edges$Target), unique(edges$TF))
allNodes <- intersect(allNodes, rownames(matBin))  # Filter to genes present in expression matrix
targets <- allNodes
total <- length(targets)

if (total == 0) {
  stop("No target genes found in edge list")
}

message("Processing ", total, " genes (targets + TFs for complete network coverage)...")

# Initialize results storage
boolRules <- list()

# Main processing loop
pb <- txtProgressBar(min = 0, max = length(targets), style = 3)

for (i in seq_along(targets)) {
  target <- targets[i]
  setTxtProgressBar(pb, i)
  
  tryCatch({
    # Select optimal regulators using SCENIC metadata
    regulator_info <- selectOptimalRegulatorsForTarget(target, edges, matBin, config$boolMaxRegulators)
    
    if (regulator_info$n_selected == 0) {
      # No regulators available - use intelligent fallback
      boolRules[[target]] <- createIntelligentFallback(target, regulator_info, matBin, config)
    } else {
      # Multi-method inference: templates + empirical
      best_rule <- inferBestBooleanRule(target, regulator_info, matBin, cellOrder, config)
      
      if (!is.null(best_rule) && best_rule$score >= config$boolMinScore) {
        boolRules[[target]] <- best_rule
      } else {
        # Fallback with regulator information
        boolRules[[target]] <- createIntelligentFallback(target, regulator_info, matBin, config)
      }
    }
  }, error = function(e) {
    warning("Failed to process ", target, ": ", e$message)
    boolRules[[target]] <- createFallbackRule(target, matBin, config)
  })
}

close(pb)

message("\nPrimary rule inference complete:")
message("  - Rules synthesized: ", length(boolRules), " target genes")

# --- STAGE 5: Add Self-Activation Rules for Missing Genes ---
message("=== STAGE 5: Adding self-activation rules for genes without rules ===")

# Find all genes mentioned as regulators but lacking rules
ruleGenes         <- names(boolRules)
genesInRules      <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
allMentionedGenes <- unique(c(ruleGenes, genesInRules))
genesNeedingRules <- setdiff(intersect(allMentionedGenes, rownames(matBin)), ruleGenes)

message("Found ", length(genesNeedingRules), " genes without rules")

if (length(genesNeedingRules) > 0) {
  if (length(genesNeedingRules) <= 10) {
    message("Genes: ", paste(genesNeedingRules, collapse = ", "))
  } else {
    message("Examples: ", paste(head(genesNeedingRules, 5), collapse = ", "), " ...")
  }
  
  for (gene in genesNeedingRules) {
    boolRules[[gene]] <- createIntelligentFallback(gene, list(n_selected = 0), matBin, config)
  }
}

message("Total rules after adding missing genes: ", length(boolRules))

# --- STAGE 6: Quality Control and Validation ---
message("=== STAGE 6: Quality control and BoolNet validation ===")

# Comprehensive rule validation
validationResults <- performFinalValidation(boolRules, matBin, edges)

message("Validation results:")
message("  - Total rules: ", validationResults$total_rules)
message("  - Valid BoolNet syntax: ", validationResults$valid_syntax)
message("  - Mean rule quality: ", round(validationResults$mean_quality, 3))
message("  - High quality rules (≥0.75): ", validationResults$high_quality, 
        " (", round(100 * validationResults$high_quality / validationResults$total_rules, 1), "%)")

# Method breakdown
method_counts <- table(sapply(boolRules, function(x) x$method %||% "unknown"))
message("Rule methods:")
for (method in names(method_counts)) {
  message("  - ", method, ": ", method_counts[method])
}

# --- STAGE 7: Save Results and Generate Reports ---
if (config$saveResults) {
  message("=== STAGE 7: Saving results and generating reports ===")
  
  # Save main results
  saveObject(boolRules, ptPaths$booleanRules, config, "Boolean rules")
  saveObject(geneMap, ptPaths$geneMap, config, "gene mapping")
  
  # Generate comprehensive output tables
  message("Generating output tables...")
  
  # Detailed rule analysis table
  ruleAnalysis <- data.frame(
    gene = names(boolRules),
    rule_string = sapply(boolRules, function(x) x$rule %||% ""),
    score = sapply(boolRules, function(x) x$score %||% 0),
    method = sapply(boolRules, function(x) x$method %||% "unknown"),
    n_regulators = sapply(boolRules, function(x) x$n_regulators %||% 0),
    k_used = sapply(boolRules, function(x) x$k_used %||% 0),
    base_accuracy = sapply(boolRules, function(x) {
      if (!is.null(x$score_breakdown)) x$score_breakdown$base_accuracy else x$score %||% 0
    }),
    sign_consistency = sapply(boolRules, function(x) {
      if (!is.null(x$score_breakdown)) x$score_breakdown$sign_consistency else 0
    }),
    scenic_confidence = sapply(boolRules, function(x) {
      if (!is.null(x$score_breakdown)) x$score_breakdown$scenic_confidence else 0
    }),
    has_and = grepl("&", sapply(boolRules, function(x) x$rule %||% "")),
    has_or = grepl("\\|", sapply(boolRules, function(x) x$rule %||% "")),
    has_not = grepl("!", sapply(boolRules, function(x) x$rule %||% "")),
    stringsAsFactors = FALSE
  )
  
  write.table(ruleAnalysis, ptPaths$booleanRulesAnalysis, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Flat regulator table for compatibility
  rulesFlat <- data.frame(
    gene = rep(names(boolRules), sapply(boolRules, function(x) max(1, length(x$regulators)))),
    regulator = unlist(lapply(boolRules, function(x) {
      if (length(x$regulators) == 0) return(NA)
      return(x$regulators)
    })),
    rule_quality = rep(sapply(boolRules, function(x) x$score %||% 0), 
                       sapply(boolRules, function(x) max(1, length(x$regulators)))),
    rule_string = rep(sapply(boolRules, function(x) x$rule %||% ""), 
                      sapply(boolRules, function(x) max(1, length(x$regulators)))),
    method = rep(sapply(boolRules, function(x) x$method %||% "unknown"),
                 sapply(boolRules, function(x) max(1, length(x$regulators)))),
    stringsAsFactors = FALSE
  )
  
  write.table(rulesFlat, ptPaths$booleanRulesFlat, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Generate report with SCENIC integration
  tryCatch({
    generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
    message("SCENIC-integrated report generated successfully!")
  }, error = function(e) {
    warning("Could not generate report: ", e$message)
  })
}

# --- FINAL SUMMARY ---
message("PRODUCTION BOOLEAN RULE INFERENCE COMPLETE")
message("SUCCESS: Generated ", length(boolRules), " Boolean rules using SCENIC-integrated approach")
message("Rule quality summary:")
message("  - Mean score: ", round(mean(sapply(boolRules, function(x) x$score %||% 0)), 3))
message("  - High quality (≥0.75): ", sum(sapply(boolRules, function(x) x$score %||% 0) >= 0.75), "/", length(boolRules))
message("  - Excellent quality (≥0.9): ", sum(sapply(boolRules, function(x) x$score %||% 0) >= 0.9), "/", length(boolRules))
message("Done!")