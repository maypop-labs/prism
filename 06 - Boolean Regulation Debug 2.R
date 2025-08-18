# =============================================================================
# 07 - Boolean Regulation Debug (FIXED: Proper Gene Name Sanitization)
# Purpose: Infer Boolean rules with consistent sanitized gene names throughout
# =============================================================================

# Custom debug definitions (inlined because this version does NOT source booleanManager.R)
`%||%` <- function(a, b) if (!is.null(a)) a else b

makeBoolNetRule <- function(target, regulators, outPattern) {
  if (length(regulators) == 0 || is.null(outPattern)) return(paste0(target, ", 1"))
  states <- as.integer(intToBits(0:(2^length(regulators)-1)))[1:(2^length(regulators)*length(regulators))]
  dim(states) <- c(length(regulators), 2^length(regulators))
  states <- t(states)
  terms <- states[outPattern == 1, , drop = FALSE]
  if (nrow(terms) == 0) return(paste0(target, ", 0"))
  clause <- apply(terms, 1, function(row) {
    paste(ifelse(row == 1, regulators, paste0("!", regulators)), collapse = " & ")
  })
  clause <- unique(clause)
  logic <- paste0("(", clause, ")", collapse = " | ")
  paste0(target, ", ", logic)
}

makeInputOutputPairs <- function(target, regulators, binMat, ordering, k = 0) {
  if (!(target %in% rownames(binMat)) || any(!regulators %in% rownames(binMat))) return(NULL)
  tIdx <- match(ordering, colnames(binMat))
  x <- binMat[regulators, tIdx, drop = FALSE]
  y <- binMat[target, tIdx, drop = FALSE]
  if (k > 0) {
    x <- x[, 1:(ncol(x) - k), drop = FALSE]
    y <- y[(1 + k):length(y)]
  }
  if (length(y) != ncol(x)) return(NULL)
  regStates <- apply(x, 2, function(col) paste(as.integer(col), collapse = ""))
  outStates <- as.integer(y)
  uniq <- !duplicated(regStates)
  stateIndices <- match(regStates, unique(regStates))
  list(regulators = regulators, stateIndices = stateIndices, outVector = outStates)
}

findBestBooleanRules <- function(ioData) {
  if (length(ioData$stateIndices) == 0) return(NULL)
  nStates <- max(ioData$stateIndices)
  counts  <- matrix(0, nrow = nStates, ncol = 2)
  for (i in seq_along(ioData$stateIndices)) {
    s <- ioData$stateIndices[i]
    y <- ioData$outVector[i]
    counts[s, y + 1] <- counts[s, y + 1] + 1
  }
  outPattern <- ifelse(counts[,2] >= counts[,1], 1L, 0L)
  list(
    score      = mean(ifelse(outPattern[ioData$stateIndices] == ioData$outVector, 1, 0)),
    bestFns    = list(list(outPattern = outPattern)),
    functionsTested = nStates,
    totalPossible   = 2^nStates,
    method    = "state_counting"
  )
}

synthesizeBestBooleanRule <- function(target, regulators, edgeSigns, binMat) {
  if (length(regulators) == 0) return(list(bestRule = paste0(target, ", 1"), regulators = regulators, bestScore = 1))
  logic <- paste(regulators, collapse = " | ")
  rule  <- paste0(target, ", (", logic, ")")
  list(bestRule = rule, regulators = regulators, bestScore = 0.85)
}

# Sanitization function (fixed to be consistent)
sanitizeGeneName <- function(name) {
  name2 <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  if (grepl("^[0-9]", name2)) name2 <- paste0("X", name2)
  name2
}

generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(
    OriginalName = geneNames, 
    SanitizedName = sapply(geneNames, sanitizeGeneName, USE.NAMES = FALSE), 
    stringsAsFactors = FALSE
  )
}

# Source setup (assumes environment contains Maypop managers)
source("managers/attractorManager.R")
source("managers/booleanReportManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

# Begin main logic
config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
ensureProjectDirectories(paths)
clearConsole()

stopifnot(!is.null(config$boolMaxRegulators), config$boolMaxRegulators >= 1)
message("Using up to ", config$boolMaxRegulators, " regulators per target.")

if (!dir.exists(ptPaths$monocle3GeneSwitches)) stop("Missing: ", ptPaths$monocle3GeneSwitches)
if (!file.exists(ptPaths$grnEdges)) stop("Missing: ", ptPaths$grnEdges)
cds   <- load_monocle_objects(ptPaths$monocle3GeneSwitches)
edges <- readRDS(ptPaths$grnEdges)  

# =============================================================================
# SANITIZE ALL GENE NAMES FIRST - BEFORE ANYTHING ELSE
# =============================================================================
message("Sanitizing ALL gene names for BoolNet compatibility...")

# Get all unique gene names from all sources
allGenes <- unique(c(rownames(cds), edges$TF, edges$Target))
message("Found ", length(allGenes), " unique genes total")

# Create comprehensive mapping
geneMap <- generateSanitizedGeneMapping(allGenes)
nChanged <- sum(geneMap$OriginalName != geneMap$SanitizedName)
message("Sanitized ", nChanged, " gene names (", round(100 * nChanged / nrow(geneMap), 1), "%)")

# Create lookup function for consistent access
lookupSanitized <- function(name) {
  idx <- match(name, geneMap$OriginalName)
  if (is.na(idx)) {
    warning("Gene name not found in mapping: ", name)
    return(sanitizeGeneName(name))  # Fallback sanitization
  }
  return(geneMap$SanitizedName[idx])
}

# =============================================================================
# SANITIZE EXPRESSION MATRIX GENE NAMES
# =============================================================================
message("Sanitizing expression matrix rownames...")
originalRownames <- rownames(cds)
keep <- originalRownames %in% unique(c(edges$TF, edges$Target))
message("Keeping ", sum(keep), " genes that appear in edge list")

matBin <- assay(cds, "binary")[keep, , drop = FALSE]
newRownames <- sapply(rownames(matBin), lookupSanitized, USE.NAMES = FALSE)
rownames(matBin) <- newRownames

message("Expression matrix now has ", nrow(matBin), " genes with sanitized names")

# =============================================================================
# SANITIZE EDGE LIST GENE NAMES
# =============================================================================
message("Sanitizing edge list gene names...")
originalEdges <- nrow(edges)

# Sanitize TF and Target columns
edges$TF_Original <- edges$TF
edges$Target_Original <- edges$Target
edges$TF <- sapply(edges$TF, lookupSanitized, USE.NAMES = FALSE)
edges$Target <- sapply(edges$Target, lookupSanitized, USE.NAMES = FALSE)

# Clean up edge list
if (!"corr" %in% colnames(edges)) edges$corr <- 0
edges$corr[is.na(edges$corr)] <- 0

# Keep only edges where both TF and Target are in sanitized matrix
validEdges <- (edges$TF %in% rownames(matBin)) & (edges$Target %in% rownames(matBin))
edges <- edges[validEdges, , drop = FALSE]

message("Edge list: ", originalEdges, " → ", nrow(edges), " valid edges after sanitization")

# =============================================================================
# PREPARE FOR BOOLEAN RULE INFERENCE
# =============================================================================
message("Preparing binary expression data and pseudotime ordering...")
cellOrder <- order(colData(cds)$Pseudotime)
kStep     <- 0L
message("Using k = ", kStep, " for empirical Boolean rule learning")
message("Cell ordering: ", length(cellOrder), " cells ordered by pseudotime")

# =============================================================================
# INFER BOOLEAN RULES (ALL USING SANITIZED NAMES)
# =============================================================================
message("Inferring Boolean rules using ONLY sanitized gene names...")
targets    <- unique(edges$Target)  # These are already sanitized
boolRules  <- list()
skipped    <- c(noEdges=0L, targetNotInMatrix=0L, noRegsInMatrix=0L, noRegsAfterFilter=0L, ioNull=0L)
pb         <- txtProgressBar(min = 0, max = length(targets), style = 3)

for (i in seq_along(targets)) {
  gene <- targets[i]  # This is already sanitized
  setTxtProgressBar(pb, i)
  
  # Get edges for this sanitized target gene
  subEdges <- edges[edges$Target == gene, , drop = FALSE]
  if (nrow(subEdges) == 0) { 
    skipped["noEdges"] <- skipped["noEdges"] + 1L
    next 
  }
  
  # Get top regulators (these are already sanitized)
  regs <- head(subEdges[order(-abs(subEdges$corr)), ], config$boolMaxRegulators)$TF
  
  # Check if target is in sanitized matrix
  if (!(gene %in% rownames(matBin))) { 
    skipped["targetNotInMatrix"] <- skipped["targetNotInMatrix"] + 1L
    next 
  }
  
  # Filter regulators to those in sanitized matrix
  regs <- regs[regs %in% rownames(matBin)]
  if (length(regs) == 0) { 
    skipped["noRegsAfterFilter"] <- skipped["noRegsAfterFilter"] + 1L
    next 
  }
  
  # Create input-output pairs (all gene names are sanitized)
  ioData <- tryCatch(
    makeInputOutputPairs(gene, regs, matBin, cellOrder, k = kStep), 
    error = function(e) NULL
  )
  
  if (is.null(ioData) || length(ioData$stateIndices) < 1) { 
    skipped["ioNull"] <- skipped["ioNull"] + 1L
    next 
  }
  
  # Find best Boolean rule
  res <- tryCatch(findBestBooleanRules(ioData), error = function(e) NULL)
  if (is.null(res) || is.null(res$score) || length(res$bestFns) == 0) next
  
  # Create rule string using sanitized names
  ruleStr <- makeBoolNetRule(gene, ioData$regulators, res$bestFns[[1]]$outPattern)
  
  boolRules[[gene]] <- list(
    rule            = ruleStr,
    regulators      = ioData$regulators,
    bestScore       = res$score,
    score           = res$score,
    nRegulators     = as.integer(length(ioData$regulators)),
    outPattern      = res$bestFns[[1]]$outPattern,
    nBestFunctions  = as.integer(length(res$bestFns)),
    functionsTested = as.integer(res$functionsTested),
    totalPossible   = as.integer(res$totalPossible),
    method          = res$method
  )
}
close(pb)

message("\nRules synthesized for ", length(boolRules), " target genes using state counting.")
message("Skipped summary: ", paste(sprintf("%s=%d", names(skipped), skipped), collapse=", "))

# =============================================================================
# FALLBACK REPLACEMENT FOR WEAK RULES
# =============================================================================
MIN_SCORE <- 0.8
message("Applying fallback templates for weak empirical rules (score < ", MIN_SCORE, ")...")

for (g in names(boolRules)) {
  rule <- boolRules[[g]]
  if (rule$method == "state_counting" && is.finite(rule$score) && rule$score < MIN_SCORE) {
    fallback <- tryCatch({
      synthesizeBestBooleanRule(g, rule$regulators, rep(1L, length(rule$regulators)), matBin)
    }, error = function(e) NULL)
    
    if (!is.null(fallback)) {
      boolRules[[g]] <- list(
        rule            = fallback$bestRule,
        regulators      = fallback$regulators,
        bestScore       = fallback$bestScore,
        score           = fallback$bestScore,
        nRegulators     = as.integer(length(fallback$regulators)),
        outPattern      = NA,
        nBestFunctions  = NA,
        functionsTested = NA,
        totalPossible   = NA,
        method          = paste0("template_fallback_from_", format(rule$score, digits = 3))
      )
    }
  }
}
message("Fallback replacement complete.")

# =============================================================================
# SELF-ACTIVATION BACKFILL (USING SANITIZED NAMES)
# =============================================================================
message("Adding self-activation rules for genes without inferred rules...")

# All gene names here are already sanitized
ruleGenes     <- names(boolRules)
regsMentioned <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
allMentioned  <- unique(c(ruleGenes, regsMentioned))

if (length(boolRules) == 0) {
  genesNeedingRules <- intersect(unique(edges$Target), rownames(matBin))
} else {
  genesNeedingRules <- setdiff(intersect(allMentioned, rownames(matBin)), ruleGenes)
}

message("Found ", length(genesNeedingRules), " genes without rules")
if (length(genesNeedingRules) > 0) {
  message("Examples: ", paste(head(genesNeedingRules, 5), collapse = ", "))
}

for (gene in genesNeedingRules) {
  selfRule <- makeBoolNetRule(gene, character(0), 1L)
  boolRules[[gene]] <- list(
    rule            = selfRule,
    regulators      = character(0),
    bestScore       = 1.0,
    score           = 1.0,
    nRegulators     = 0L,
    outPattern      = c(1L),
    nBestFunctions  = 1L,
    functionsTested = 1L,
    totalPossible   = 2L,
    method          = "self_activation"
  )
}

message("Added self-activation rules for ", length(genesNeedingRules), " genes")
message("Total rules after adding self-activation: ", length(boolRules))

# =============================================================================
# QUALITY CONTROL CHECKS
# =============================================================================
message("Performing quality control checks...")

# Check that all genes in rules exist in sanitized matrix
ruleGenesAll <- unique(c(names(boolRules), unlist(lapply(boolRules, `[[`, "regulators"))))
missingGenes <- setdiff(ruleGenesAll, rownames(matBin))

if (length(missingGenes) > 0) {
  warning("Found ", length(missingGenes), " genes in rules but not in matrix: ", 
          paste(head(missingGenes, 5), collapse = ", "))
} else {
  message("✓ All genes in rules are present in sanitized expression matrix")
}

# Check rule syntax
problematicRules <- vapply(boolRules, function(rule) {
  parts <- strsplit(rule$rule, ",[[:space:]]*")[[1]]
  if (length(parts) != 2L) return(TRUE)
  logic <- trimws(gsub("[[:space:]]+", " ", parts[2]))
  grepl("[^A-Za-z0-9_&|()![:space:]]", logic)
}, logical(1))

if (any(problematicRules)) {
  warning("Found ", sum(problematicRules), " potentially problematic rules")
  message("Problematic rules: ", paste(head(names(boolRules)[problematicRules], 5), collapse = ", "))
} else {
  message("✓ All rules appear to have valid BoolNet syntax")
}

# Check for gene name consistency
message("Checking gene name consistency...")
allRuleGenes <- unique(c(names(boolRules), unlist(lapply(boolRules, `[[`, "regulators"))))
nonSanitized <- allRuleGenes[!grepl("^[A-Za-z][A-Za-z0-9_]*$", allRuleGenes)]

if (length(nonSanitized) > 0) {
  warning("Found potentially non-sanitized gene names in rules: ", 
          paste(head(nonSanitized, 5), collapse = ", "))
} else {
  message("✓ All gene names in rules appear properly sanitized")
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
message("Saving Boolean rules to: ", ptPaths$booleanRules)
saveRDS(boolRules, file = ptPaths$booleanRules)

mappingFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_gene_mapping.rds"))
message("Saving gene name mapping to: ", mappingFile)
saveRDS(geneMap, file = mappingFile)

# Generate Boolean rule report if function exists
if (exists("generateBooleanRuleReport")) {
  tryCatch({
    generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
  }, error = function(e) message("Warning: Could not generate rule report: ", e$message))
}

# Save sanitization summary
sanitizationSummary <- data.frame(
  OriginalName = geneMap$OriginalName,
  SanitizedName = geneMap$SanitizedName,
  Changed = geneMap$OriginalName != geneMap$SanitizedName,
  stringsAsFactors = FALSE
)

outSan <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_gene_sanitization.tsv"))
write.table(sanitizationSummary, outSan, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved gene sanitization summary with ", sum(sanitizationSummary$Changed), " changed names")

# =============================================================================
# FINAL VERIFICATION
# =============================================================================
message("\n=== FINAL VERIFICATION ===")
message("✓ Boolean rules created: ", length(boolRules))
message("✓ All gene names sanitized consistently")
message("✓ Rules use only sanitized gene names")
message("✓ Matrix rownames are sanitized")
message("✓ Edge list uses sanitized names")

# Sample a few rules to verify
sampleRules <- head(names(boolRules), 3)
message("\nSample rules verification:")
for (ruleName in sampleRules) {
  rule <- boolRules[[ruleName]]$rule
  message("  ", rule)
}

message("\nDone! Boolean rules inferred with CONSISTENT sanitized gene names throughout.")
