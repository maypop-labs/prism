# =============================================================================
# 06 - Boolean Regulation (Production Version)
# Purpose: Generate production-quality Boolean rules using k=0 + SCENIC data
# =============================================================================

# Core Philosophy: 
# - k = 0 only (same pseudotime, maximum statistical power)
# - Leverage all SCENIC/GENIE3 information (regType, motif confidence, correlation)
# - Multi-method approach (empirical + template + fallback)
# - Bulletproof error handling and resumable execution
# - BoolNet-ready output with comprehensive validation

source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/booleanReportManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

library(BoolNet)
library(dplyr)
library(Matrix)

# --- CONFIGURATION ---
config <- initializeScript()
config$boolMaxRegulators <- 3
config$boolMinScore <- 0.05
config$boolMinStates <- 10
config$kValue <- 0  # Fixed at 0 for production
config$checkpointInterval <- 25
config$weightsConfig <- list(
  correlation = 1.0,
  motif_confidence = 0.3,
  scenic_importance = 0.2,
  prior_knowledge = 0.1
)

pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
ensureProjectDirectories(paths)
clearConsole()

# --- STAGE 1: Robust Data Loading and Validation ---
message("=== STAGE 1: Loading and validating data ===")

if (!dir.exists(ptPaths$monocle3GeneSwitches)) {
  stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
}
if (!file.exists(ptPaths$grnEdges)) {
  stop("Edges RDS file not found: ", ptPaths$grnEdges)
}

cds   <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)
edges <- readRDS(ptPaths$grnEdges)

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
if (!"motifConfidence" %in% colnames(edges)) {
  edges$motifConfidence <- NA
}
if (!"NES" %in% colnames(edges)) {
  edges$NES <- NA
}
if (!"hasMotif" %in% colnames(edges)) {
  edges$hasMotif <- NA
}

message("Data validation complete:")
message("  - CDS object: ", nrow(cds), " genes, ", ncol(cds), " cells")
message("  - Edge list: ", nrow(edges), " regulatory relationships")
message("  - SCENIC columns: ", paste(colnames(edges), collapse = ", "))

# --- STAGE 2: Comprehensive Gene Name Sanitization ---
message("=== STAGE 2: Sanitizing gene names for BoolNet compatibility ===")

# Generate comprehensive gene name mapping
allGenes <- unique(c(rownames(cds), edges$TF, edges$Target))
geneMap  <- generateSanitizedGeneMapping(allGenes)
geneLookup <- setNames(geneMap$SanitizedName, geneMap$OriginalName)

message("Gene name mapping:")
message("  - Total unique genes: ", length(allGenes))
message("  - Genes requiring sanitization: ", sum(geneMap$OriginalName != geneMap$SanitizedName))

# Apply sanitization to expression matrix
matBin <- assay(cds, "binary")
genesInEdges <- unique(c(edges$TF, edges$Target))
matBin <- matBin[rownames(matBin) %in% genesInEdges, ]

originalRowNames <- rownames(matBin)
sanitizedRowNames <- sapply(originalRowNames, function(g) geneLookup[[g]] %||% g)
rownames(matBin) <- sanitizedRowNames

# Apply sanitization to edge list
edges$TF <- sapply(edges$TF, function(g) geneLookup[[g]] %||% g)
edges$Target <- sapply(edges$Target, function(g) geneLookup[[g]] %||% g)

# Remove edges with genes missing from expression matrix
validTFs <- edges$TF %in% rownames(matBin)
validTargets <- edges$Target %in% rownames(matBin)
validEdges <- validTFs & validTargets

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
kStep <- 0  # Fixed at k=0 for production

message("Configuration:")
message("  - Using k = ", kStep, " (same pseudotime analysis)")
message("  - Cell ordering: ", length(cellOrder), " cells ordered by pseudotime")
message("  - Max regulators per gene: ", config$boolMaxRegulators)
message("  - Minimum rule score: ", config$boolMinScore)

# --- STAGE 4: Enhanced Boolean Rule Inference ---
message("=== STAGE 4: SCENIC-enhanced Boolean rule inference ===")

targets <- unique(edges$Target)
total <- length(targets)

if (total == 0) {
  stop("No target genes found in edge list")
}

message("Processing ", total, " target genes...")

# Initialize results storage
boolRules <- list()

# Resume support
checkpoint_file <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_boolean_checkpoint.rds"))
completed_genes <- character(0)

if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  boolRules <- checkpoint$rules
  completed_genes <- checkpoint$completed
  targets <- setdiff(targets, completed_genes)
  message("Resuming from checkpoint: ", length(completed_genes), " genes already completed")
}

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
    
    # Checkpoint every N genes
    if (i %% config$checkpointInterval == 0) {
      completed_genes <- c(completed_genes, names(boolRules))
      saveRDS(list(rules = boolRules, completed = unique(completed_genes)), checkpoint_file)
      message("\nCheckpoint saved at gene ", i, "/", length(targets))
    }
    
  }, error = function(e) {
    warning("Failed to process ", target, ": ", e$message)
    boolRules[[target]] <- createFallbackRule(target, matBin, config)
  })
}

close(pb)

# Final checkpoint cleanup
if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
}

message("\nPrimary rule inference complete:")
message("  - Rules synthesized: ", length(boolRules), " target genes")

# --- STAGE 5: Add Self-Activation Rules for Missing Genes ---
message("=== STAGE 5: Adding self-activation rules for genes without rules ===")

# Find all genes mentioned as regulators but lacking rules
ruleGenes <- names(boolRules)
genesInRules <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
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
  saveRDS(boolRules, file = ptPaths$booleanRules)
  message("Boolean rules saved to: ", ptPaths$booleanRules)
  
  # Save gene name mapping
  mappingFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_gene_mapping.rds"))
  saveRDS(geneMap, file = mappingFile)
  message("Gene name mapping saved to: ", mappingFile)
  
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
  
  analysisFile <- paste0(paths$base$tsv, cellType, "_", trajectory, "_boolean_rules_analysis.tsv")
  write.table(ruleAnalysis, analysisFile, sep = "\t", row.names = FALSE, quote = FALSE)
  
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
  
  flatFile <- paste0(paths$base$tsv, cellType, "_", trajectory, "_boolean_rules.tsv")
  write.table(rulesFlat, flatFile, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Gene sanitization summary
  sanitizationSummary <- data.frame(
    OriginalName = geneMap$OriginalName,
    SanitizedName = geneMap$SanitizedName,
    Changed = geneMap$OriginalName != geneMap$SanitizedName,
    stringsAsFactors = FALSE
  )
  
  sanitFile <- paste0(paths$base$tsv, cellType, "_", trajectory, "_gene_sanitization.tsv")
  write.table(sanitizationSummary, sanitFile, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("Output tables saved:")
  message("  - Analysis table: ", basename(analysisFile))
  message("  - Flat rules table: ", basename(flatFile))
  message("  - Sanitization summary: ", basename(sanitFile))
  message("  - Genes sanitized: ", sum(sanitizationSummary$Changed))
  
  # Generate enhanced report with SCENIC integration
  tryCatch({
    generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
    message("Enhanced SCENIC-integrated report generated successfully!")
  }, error = function(e) {
    warning("Could not generate enhanced report: ", e$message)
  })
}

# --- FINAL SUMMARY ---
message("\n" , paste(rep("=", 60), collapse = ""))
message("PRODUCTION BOOLEAN RULE INFERENCE COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("SUCCESS: Generated ", length(boolRules), " Boolean rules using SCENIC-enhanced approach")
message("")
message("Key achievements:")
message("✓ Used k=0 (same pseudotime) for maximum statistical power")
message("✓ Integrated SCENIC metadata (regType, motif confidence, NES scores)")
message("✓ Multi-method approach: templates + empirical + intelligent fallbacks")
message("✓ Comprehensive quality control and BoolNet validation")
message("✓ Production-grade error handling and resumable execution")
message("")
message("Rule quality summary:")
message("  - Mean score: ", round(mean(sapply(boolRules, function(x) x$score %||% 0)), 3))
message("  - High quality (≥0.75): ", sum(sapply(boolRules, function(x) x$score %||% 0) >= 0.75), 
        "/", length(boolRules))
message("  - Excellent quality (≥0.9): ", sum(sapply(boolRules, function(x) x$score %||% 0) >= 0.9), 
        "/", length(boolRules))
message("")
message("Ready for attractor analysis and perturbation experiments!")
message(paste(rep("=", 60), collapse = ""))
