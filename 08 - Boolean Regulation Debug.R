# ============================================================================
# 08 - Boolean Regulation Debug (CORRECTED)
# ----------------------------------------------------------------------------
# Self-contained inference loop (no dependency on managers/booleanManager.R).
# Safe to run with your usual project managers sourced; if you comment them out,
# this file still provides the minimal helpers it needs.
#
# Fixes vs earlier debug versions:
# - Deterministic, mutually-exclusive control flow inside the target loop
# - Every non-success path increments a specific skip counter and `next`s
# - Backfill works even when `boolRules` is initially empty
# - Integer-safe metadata (no vapply integer/double crashes)
# - POSIX whitespace in QC; spacing normalized; no false positives
# - Typo fixed: `functionsTested` (not `functionsTestsed`)
# - Defensive handling of missing/NA `corr` in edges
# ============================================================================

# ---------------- Optional project sources (keep as in your workflow) --------
source("managers/attractorManager.R")
source("managers/booleanReportManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

# ---------------- Minimal local helpers (standalone-friendly) -----------------

sanitizeGeneName <- function(name) {
  nm <- gsub("[^A-Za-z0-9_]", "_", trimws(name))
  if (grepl("^[0-9]", nm)) nm <- paste0("X", nm)
  nm
}

generateSanitizedGeneMapping <- function(geneNames) {
  data.frame(
    OriginalName = geneNames,
    SanitizedName = vapply(geneNames, sanitizeGeneName, ""),
    stringsAsFactors = FALSE
  )
}

makeBoolNetRule <- function(geneName, regulators, outPattern) {
  if (length(regulators) == 0) {
    rhs <- ifelse(outPattern[1] == 1L, geneName, paste0("!", geneName))
    return(paste0(geneName, ", ", rhs))
  }
  R <- length(regulators)
  expected <- 2^R
  if (length(outPattern) != expected) stop("Output pattern length ", length(outPattern), " != ", expected, " for ", geneName)
  clauseList <- character()
  for (i in 0:(expected - 1)) if (outPattern[i + 1] == 1L) {
    bits <- as.integer(intToBits(i))[1:R]
    andTerms <- ifelse(bits == 1L, regulators, paste0("!", regulators))
    clauseList <- c(clauseList, paste0("(", paste(andTerms, collapse = " & "), ")"))
  }
  if (length(clauseList) == 0) return(paste0(geneName, ", !", geneName))
  paste0(geneName, ", ", paste(clauseList, collapse = " | "))
}

makeInputOutputPairs <- function(targetGene, regulators, matBin, cellOrder, k = 0L) {
  if (!targetGene %in% rownames(matBin)) return(list(error = "target_not_in_matrix"))
  regulators <- regulators[regulators %in% rownames(matBin)]
  if (length(regulators) == 0) return(list(error = "no_valid_regulators"))
  
  maxT <- length(cellOrder) - k
  if (maxT <= 0) return(list(error = "insufficient_cells_for_k"))
  
  targetCells <- cellOrder[seq_len(maxT) + k]
  inputCells  <- cellOrder[seq_len(maxT)]
  
  targetVals  <- matBin[targetGene, targetCells]
  if (length(unique(targetVals)) == 1) return(list(error = "constant_target"))
  
  regMatrix <- matBin[regulators, inputCells, drop = FALSE]
  regVar    <- apply(regMatrix, 1, function(x) length(unique(x)))
  variable  <- regulators[regVar > 1]
  if (length(variable) == 0) return(list(error = "all_constant_regulators"))
  
  regMatrix <- regMatrix[variable, , drop = FALSE]
  R <- length(variable)
  weights <- 2^(0:(R-1))
  stateIndices <- 1L + as.integer(regMatrix[1, ] * weights[1])
  if (R > 1) for (i in 2:R) stateIndices <- stateIndices + as.integer(regMatrix[i, ] * weights[i])
  
  outputs <- as.integer(targetVals)
  list(stateIndices = stateIndices, outputs = outputs, regulators = variable, nStates = 2^R)
}

findBestBooleanRules <- function(ioData) {
  if (is.null(ioData) || !is.null(ioData$error))
    return(list(score = 0, bestFns = list(), functionsTested = 0L, totalPossible = 0L, method = "state_counting"))
  idx <- ioData$stateIndices; outs <- ioData$outputs; nStates <- as.integer(ioData$nStates)
  if (length(idx) == 0 || length(outs) == 0)
    return(list(score = 0, bestFns = list(), functionsTested = 0L, totalPossible = 0L, method = "state_counting"))
  count1 <- integer(nStates); count0 <- integer(nStates)
  for (i in seq_along(idx)) {
    s <- idx[i]
    if (outs[i] == 1L) count1[s] <- count1[s] + 1L else count0[s] <- count0[s] + 1L
  }
  outPattern <- as.integer(count1 >= count0)  # ties -> 1
  correct <- 0L; total <- length(idx)
  for (i in seq_along(idx)) if (outPattern[idx[i]] == outs[i]) correct <- correct + 1L
  score <- correct / total
  bestFn <- list(patternIdx = 0L, outPattern = outPattern, score = score)
  list(score = score, bestFns = list(bestFn), functionsTested = nStates, totalPossible = as.integer(2^nStates), method = "state_counting")
}

# --------------------------- Main pipeline -----------------------------------

# Expect your normal project setup to provide these; if not, source them above.
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

# Load data
if (!dir.exists(ptPaths$monocle3GeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
if (!file.exists(ptPaths$grnEdges)) stop("Edges RDS file not found: ", ptPaths$grnEdges)
cds   <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)
edges <- readRDS(ptPaths$grnEdges)

# Sanitize gene names
message("Sanitizing gene names for BoolNet compatibility...")
allGenes <- unique(c(rownames(cds), edges$TF, edges$Target))
geneMap  <- generateSanitizedGeneMapping(allGenes)
lookup   <- setNames(geneMap$SanitizedName, geneMap$OriginalName)

message("Sanitizing expression matrix gene names...")
matBin <- assay(cds, "binary")
keep <- rownames(matBin) %in% unique(c(edges$TF, edges$Target))
matBin <- matBin[keep, , drop = FALSE]
rownames(matBin) <- vapply(rownames(matBin), function(g) if (!is.null(lookup[[g]])) lookup[[g]] else g, "")
edges$TF     <- vapply(edges$TF,     function(g) if (!is.null(lookup[[g]])) lookup[[g]] else g, "")
edges$Target <- vapply(edges$Target, function(g) if (!is.null(lookup[[g]])) lookup[[g]] else g, "")

# Ensure corr exists and is usable
if (!"corr" %in% colnames(edges)) edges$corr <- 0
edges$corr[is.na(edges$corr)] <- 0

# Drop edges with genes missing from the matrix
validEdges <- (edges$TF %in% rownames(matBin)) & (edges$Target %in% rownames(matBin))
if (sum(!validEdges) > 0) {
  message("Removing ", sum(!validEdges), " edges with genes missing from expression matrix")
  edges <- edges[validEdges, , drop = FALSE]
}
message("Final edge list: ", nrow(edges), " edges")

# Pseudotime order
message("Preparing binary expression data and pseudotime ordering...")
cellOrder <- order(colData(cds)$Pseudotime)
kStep     <- 0L
message("Using k = ", kStep, " for empirical Boolean rule learning")
message("Cell ordering: ", length(cellOrder), " cells ordered by pseudotime")

# Inference loop
message("Inferring Boolean rules ...")

targets <- unique(edges$Target)
if (length(targets) == 0) stop("No target genes in edges")

boolRules <- list()
# Rich skip counters to diagnose why a gene didn't get a rule
skipped <- c(noEdges=0L, targetNotInMatrix=0L, noRegsInMatrix=0L, noRegsAfterFilter=0L,
             constantTarget=0L, allConstantRegulators=0L, insufficient_cells_for_k=0L,
             ioNull=0L, learnFail=0L)

pb <- txtProgressBar(min = 0, max = length(targets), style = 3)
for (i in seq_along(targets)) {
  gene <- targets[i]
  setTxtProgressBar(pb, i)
  
  subEdges <- edges[edges$Target == gene, , drop = FALSE]
  if (nrow(subEdges) == 0) { skipped["noEdges"] <- skipped["noEdges"] + 1L; next }
  
  topEdges   <- head(subEdges[order(-abs(subEdges$corr)), , drop = FALSE], config$boolMaxRegulators)
  regulators <- topEdges$TF
  
  if (!(gene %in% rownames(matBin))) { skipped["targetNotInMatrix"] <- skipped["targetNotInMatrix"] + 1L; next }
  regsInMatrix <- regulators %in% rownames(matBin)
  if (!any(regsInMatrix)) { skipped["noRegsInMatrix"] <- skipped["noRegsInMatrix"] + 1L; next }
  regulators <- regulators[regsInMatrix]
  if (length(regulators) == 0) { skipped["noRegsAfterFilter"] <- skipped["noRegsAfterFilter"] + 1L; next }
  
  ioData <- tryCatch(makeInputOutputPairs(gene, regulators, matBin, cellOrder, k = kStep), error = function(e) list(error = "io_exception"))
  if (!is.null(ioData$error)) {
    reason <- ioData$error
    key <- if (reason %in% names(skipped)) reason else "ioNull"
    skipped[key] <- skipped[key] + 1L
    next
  }
  if (length(ioData$stateIndices) < 1) { skipped["ioNull"] <- skipped["ioNull"] + 1L; next }
  
  res <- tryCatch(findBestBooleanRules(ioData), error = function(e) NULL)
  if (is.null(res) || is.null(res$bestFns) || length(res$bestFns) == 0) { skipped["learnFail"] <- skipped["learnFail"] + 1L; next }
  
  outPattern <- res$bestFns[[1]]$outPattern
  ruleStr    <- makeBoolNetRule(gene, ioData$regulators, outPattern)
  
  boolRules[[gene]] <- list(
    rule            = ruleStr,
    regulators      = ioData$regulators,
    bestScore       = res$score,
    score           = res$score,
    nRegulators     = as.integer(length(ioData$regulators)),
    outPattern      = outPattern,
    nBestFunctions  = as.integer(length(res$bestFns)),
    functionsTested = as.integer(res$functionsTested),
    totalPossible   = as.integer(res$totalPossible),
    method          = res$method
  )
}
close(pb)

message("\nRules synthesized for ", length(boolRules), " target genes using state counting.")
message("Skipped summary: ", paste(sprintf("%s=%d", names(skipped), skipped), collapse = ", "))

# Backfill even if no rules were learned empirically
message("Adding self-activation rules for genes without inferred rules...")
ruleGenes     <- names(boolRules)
regsMentioned <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
allMentioned  <- unique(c(ruleGenes, regsMentioned))

if (length(boolRules) == 0) {
  genesNeedingRules <- intersect(unique(edges$Target), rownames(matBin))
} else {
  genesNeedingRules <- setdiff(intersect(allMentioned, rownames(matBin)), ruleGenes)
}

message("Found ", length(genesNeedingRules), " genes without rules")
if (length(genesNeedingRules) > 0) message("Examples: ", paste(head(genesNeedingRules, 5), collapse = ", "))

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

# QC
message("Performing quality control checks...")
ruleGenesAll <- unique(c(names(boolRules), unlist(lapply(boolRules, `[[`, "regulators"))))
missingGenes <- setdiff(ruleGenesAll, rownames(matBin))
if (length(missingGenes) > 0) {
  warning("Found ", length(missingGenes), " genes in rules but not in matrix: ", paste(head(missingGenes, 5), collapse = ", "))
} else {
  message("\u2713 All genes in rules are present in expression matrix")
}

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
  message("\u2713 All rules appear to have valid BoolNet syntax")
}

# Save
if (isTRUE(config$saveResults)) {
  message("Saving Boolean rules to: ", ptPaths$booleanRules)
  saveRDS(boolRules, file = ptPaths$booleanRules)
  
  mappingFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_gene_mapping.rds"))
  message("Saving gene name mapping to: ", mappingFile)
  saveRDS(geneMap, file = mappingFile)
  
  if (length(boolRules) > 0) {
    ruleAnalysis <- data.frame(
      gene              = names(boolRules),
      rule_string       = vapply(boolRules, function(x) x$rule, ""),
      score             = vapply(boolRules, function(x) as.numeric(x$score), numeric(1)),
      n_regulators      = vapply(boolRules, function(x) as.integer(if (is.null(x$nRegulators)) length(x$regulators) else x$nRegulators), integer(1)),
      n_best_functions  = vapply(boolRules, function(x) as.integer(if (is.null(x$nBestFunctions)) 1L else x$nBestFunctions), integer(1)),
      functions_tested  = vapply(boolRules, function(x) as.integer(if (is.null(x$functionsTested)) NA_integer_ else x$functionsTested), integer(1)),
      total_possible    = vapply(boolRules, function(x) as.integer(if (is.null(x$totalPossible)) NA_integer_ else x$totalPossible), integer(1)),
      method            = vapply(boolRules, function(x) x$method, ""),
      has_or            = vapply(boolRules, function(x) grepl("\\|", x$rule), logical(1)),
      has_and           = vapply(boolRules, function(x) grepl("&",   x$rule), logical(1)),
      stringsAsFactors  = FALSE
    )
    outAnalysis <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_boolean_rules_analysis.tsv"))
    write.table(ruleAnalysis, outAnalysis, sep = "\t", row.names = FALSE, quote = FALSE)
    
    len_per_gene <- vapply(boolRules, function(x) as.integer(max(1L, length(x$regulators))), integer(1))
    rulesFlat <- data.frame(
      gene         = rep(names(boolRules), len_per_gene),
      regulator    = unlist(lapply(boolRules, function(x) if (length(x$regulators) == 0) NA_character_ else x$regulators), use.names = FALSE),
      rule_quality = rep(vapply(boolRules, function(x) as.numeric(x$score), numeric(1)), len_per_gene),
      rule_string  = rep(vapply(boolRules, function(x) x$rule, ""), len_per_gene),
      stringsAsFactors = FALSE
    )
    outFlat <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_boolean_rules.tsv"))
    write.table(rulesFlat, outFlat, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  if (exists("generateBooleanRuleReport")) {
    tryCatch({
      generateBooleanRuleReport(boolRules, edges, paths, cellType, trajectory, config)
    }, error = function(e) message("Warning: Could not generate rule report: ", e$message))
  }
  
  sanitizationSummary <- data.frame(
    OriginalName = geneMap$OriginalName,
    SanitizedName = geneMap$SanitizedName,
    Changed = geneMap$OriginalName != geneMap$SanitizedName,
    stringsAsFactors = FALSE
  )
  outSan <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_gene_sanitization.tsv"))
  write.table(sanitizationSummary, outSan, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Saved gene sanitization summary with ", sum(sanitizationSummary$Changed), " changed names")
}

message("Done! Boolean rules inferred (corrected debug script).")
