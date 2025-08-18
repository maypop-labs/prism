# =============================================================================
# 08 - BoolNet Debug.R
#
# Load Boolean rules, create BoolNet object, compute attractors, save results.
# Updated to work with optimized Boolean rules from state counting approach.
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

# --- Load Boolean rules ---
message("Loading Boolean rules from: ", ptPaths$booleanRules)
boolRules <- readRDS(ptPaths$booleanRules)
message("Loaded ", length(boolRules), " Boolean rules")

# --- Analyze rule quality ---
message("Analyzing rule quality...")
empiricalRules <- boolRules[sapply(boolRules, function(x) x$method == "state_counting")]
selfActivationRules <- boolRules[sapply(boolRules, function(x) x$method == "self_activation")]
orRules <- boolRules[sapply(boolRules, function(x) grepl("\\|", x$rule))]

message("- Empirical rules (state counting): ", length(empiricalRules))
message("- Self-activation rules: ", length(selfActivationRules))
message("- Rules with OR operations: ", length(orRules))

if (length(empiricalRules) > 0) {
  scores <- sapply(empiricalRules, function(x) x$score)
  message("- Mean empirical rule score: ", round(mean(scores), 3))
  message("- Rules with score > 0.8: ", sum(scores > 0.8))
}

# --- Build BoolNet Rule Table ---
message("Building BoolNet rule table...")
ruleLines <- c("targets, factors")

for (gene in names(boolRules)) {
  ruleLines <- c(ruleLines, boolRules[[gene]]$rule)
}

message("Rule table complete: ", length(ruleLines) - 1, " rules")

# --- Build BoolNet Rule Table ---
message("Building BoolNet rule table...")
ruleLines <- c("targets, factors")

for (gene in names(boolRules)) {
  rule <- boolRules[[gene]]$rule
  # Normalize whitespace to prevent Unicode issues
  rule <- trimws(rule)
  rule <- gsub("\\s+", " ", rule)  # Collapse multiple spaces
  ruleLines <- c(ruleLines, rule)
}

message("Rule table complete: ", length(ruleLines) - 1, " rules")

# --- Create BoolNet Object (Let BoolNet validate, not us) ---
message("Creating BoolNet object...")
tempFile <- tempfile(fileext = ".txt")
writeLines(ruleLines, con = tempFile)

tryCatch({
  boolnet <- loadNetwork(tempFile)
  boolnet$type <- "synchronous"
  message("✓ BoolNet object created successfully: ", length(boolnet$genes), " genes")
  message("✓ All ", length(ruleLines) - 1, " rules were accepted by BoolNet")
}, error = function(e) {
  message("❌ BoolNet parsing failed: ", e$message)
  message("Attempting to identify problematic rules...")
  
  # Binary search to find problematic rules
  validRules <- c("targets, factors")
  problemRules <- character()
  
  for (i in 2:length(ruleLines)) {
    testRules <- c(validRules, ruleLines[i])
    tempTestFile <- tempfile(fileext = ".txt")
    writeLines(testRules, con = tempTestFile)
    
    testResult <- tryCatch({
      testNet <- loadNetwork(tempTestFile)
      unlink(tempTestFile)
      TRUE
    }, error = function(e) {
      unlink(tempTestFile)
      FALSE
    })
    
    if (testResult) {
      validRules <- testRules
    } else {
      problemRules <- c(problemRules, ruleLines[i])
      message("Problem rule found: ", ruleLines[i])
      # Show character codes to detect Unicode issues
      chars <- utf8ToInt(ruleLines[i])
      nonAscii <- chars[chars > 127]
      if (length(nonAscii) > 0) {
        message("  Contains non-ASCII characters: ", paste(nonAscii, collapse = ", "))
      }
    }
  }
  
  if (length(validRules) > 1) {
    message("Proceeding with ", length(validRules) - 1, " valid rules")
    writeLines(validRules, con = tempFile)
    boolnet <- loadNetwork(tempFile)
    boolnet$type <- "synchronous"
  } else {
    stop("No valid rules found")
  }
})

# Clean up temp file
unlink(tempFile)

# --- Compute Attractors ---
message("Computing attractors...")

# Set seed for reproducible results
set.seed(42)  # ChatGPT5's suggestion for reproducible debugging

# Use appropriate method based on network size
nGenes <- length(boolnet$genes)
if (nGenes <= 20) {
  # Small network - use exhaustive search
  message("Small network (", nGenes, " genes) - using exhaustive search")
  attractors <- getAttractors(
    boolnet,
    method = "exhaustive",
    type = "synchronous",
    returnTable = TRUE
  )
} else {
  # Larger network - use random sampling
  message("Larger network (", nGenes, " genes) - using random sampling with seed=42")
  nSamples <- min(1000, 2^min(nGenes, 15))  # Don't go crazy with samples
  attractors <- getAttractors(
    boolnet,
    method = "random",
    startStates = nSamples,
    type = "synchronous",
    returnTable = TRUE
  )
}

# --- Analyze Results ---
nAttractors <- length(attractors$attractors)
attractorSizes <- sapply(attractors$attractors, function(att) ncol(att$involvedStates))
uniqueSizes <- unique(attractorSizes)

message("ATTRACTOR ANALYSIS:")
message("- Found ", nAttractors, " attractors")
message("- Attractor sizes: ", paste(sort(uniqueSizes), collapse = ", "))

# Basin size analysis
if (nAttractors > 0) {
  basinSizes <- sapply(attractors$attractors, function(att) att$basinSize)
  totalStates <- sum(basinSizes)
  
  message("- Total states covered: ", totalStates)
  if (nGenes <= 20) {
    expectedStates <- 2^nGenes
    message("- Expected states (2^n): ", expectedStates)
    message("- Coverage: ", round(100 * totalStates / expectedStates, 1), "%")
  }
  
  # Show largest attractors
  if (nAttractors > 1) {
    sortedIdx <- order(basinSizes, decreasing = TRUE)
    message("- Largest attractor basins:")
    for (i in 1:min(5, nAttractors)) {
      idx <- sortedIdx[i]
      message("  Attractor ", idx, ": basin size ", basinSizes[idx], 
              " (", round(100 * basinSizes[idx] / totalStates, 1), "%)")
    }
  }
}

# --- Diagnostic Assessment ---
message("DIAGNOSTIC ASSESSMENT:")

if (nAttractors > 100) {
  message("❌ WARNING: Very high number of attractors (", nAttractors, ") suggests network issues")
  message("   Possible causes:")
  message("   - Rules too restrictive (many AND-only rules)")
  message("   - Poor rule quality scores")
  message("   - Insufficient OR operations for flexibility")
} else if (nAttractors > 50) {
  message("⚠️  CAUTION: High number of attractors (", nAttractors, ")")
  message("   Network may be somewhat fragmented")
} else if (nAttractors < 5) {
  message("✅ EXCELLENT: Low number of attractors (", nAttractors, ")")
  message("   Network shows strong convergent behavior")
} else {
  message("✅ GOOD: Reasonable number of attractors (", nAttractors, ")")
  message("   Network appears well-balanced")
}

# Point attractor preference check
pointAttractors <- sum(attractorSizes == 1)
cycleAttractors <- sum(attractorSizes > 1)

message("- Point attractors (steady states): ", pointAttractors)
message("- Cycle attractors: ", cycleAttractors)

if (pointAttractors > cycleAttractors) {
  message("✅ Good: Network prefers steady states over cycles")
} else {
  message("ℹ️  Info: Network has many cyclic attractors")
}

# --- Save Results ---
if (config$saveResults) {
  message("Saving results...")
  
  message("Saving BoolNet object to: ", ptPaths$boolnet)
  saveRDS(boolnet, file = ptPaths$boolnet)
  
  message("Saving attractors to: ", ptPaths$attractors)
  saveRDS(attractors, file = ptPaths$attractors)
  
  # Save diagnostic summary
  diagnosticSummary <- list(
    nGenes = nGenes,
    nRules = length(ruleLines) - 1,
    nAttractors = nAttractors,
    attractorSizes = attractorSizes,
    basinSizes = if (nAttractors > 0) sapply(attractors$attractors, function(att) att$basinSize) else integer(0),
    pointAttractors = pointAttractors,
    cycleAttractors = cycleAttractors,
    ruleQuality = if (length(empiricalRules) > 0) mean(sapply(empiricalRules, function(x) x$score)) else NA,
    timestamp = Sys.time()
  )
  
  diagnosticFile <- file.path(paths$base$rds, paste0(cellType, "_", trajectory, "_attractor_diagnostics.rds"))
  saveRDS(diagnosticSummary, file = diagnosticFile)
  message("Diagnostic summary saved to: ", diagnosticFile)
  
  # Write human-readable summary
  summaryFile <- file.path(paths$base$txt, paste0(cellType, "_", trajectory, "_attractor_summary.txt"))
  cat("ATTRACTOR ANALYSIS SUMMARY\n", file = summaryFile)
  cat("==========================\n\n", file = summaryFile, append = TRUE)
  cat("Network Size: ", nGenes, " genes\n", file = summaryFile, append = TRUE)
  cat("Rules Used: ", length(ruleLines) - 1, "\n", file = summaryFile, append = TRUE)
  cat("Attractors Found: ", nAttractors, "\n", file = summaryFile, append = TRUE)
  cat("Point Attractors: ", pointAttractors, "\n", file = summaryFile, append = TRUE)
  cat("Cycle Attractors: ", cycleAttractors, "\n", file = summaryFile, append = TRUE)
  if (length(empiricalRules) > 0) {
    cat("Mean Rule Quality: ", round(mean(sapply(empiricalRules, function(x) x$score)), 3), "\n", 
        file = summaryFile, append = TRUE)
  }
  cat("Analysis Date: ", as.character(Sys.time()), "\n", file = summaryFile, append = TRUE)
  
  message("Summary saved to: ", summaryFile)
}

message("Done! Network analysis complete.")

# Final recommendation
if (nAttractors <= 20) {
  message("\n🎯 RECOMMENDATION: Network looks good for perturbation analysis!")
} else if (nAttractors <= 50) {
  message("\n💭 RECOMMENDATION: Network may work for perturbation analysis, monitor carefully")
} else {
  message("\n⚠️  RECOMMENDATION: Consider rule quality improvements before perturbation analysis")
}