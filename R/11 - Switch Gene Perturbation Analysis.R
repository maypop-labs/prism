# =============================================================================
# 11 - Switch Gene Perturbation Analysis.R
# Purpose: Analyze whether perturbing switch genes in the direction they 
#          naturally switch affects aging scores, and vice versa
# =============================================================================

# Load required packages and initialize paths
source("managers/attractorManager.R")
source("managers/pathManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")
config <- initializeScript()
pathInfo <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths <- pathInfo$paths
cellType <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths <- getCellTypeFilePaths(paths$base, cellType)
ptPaths <- getTrajectoryFilePaths(paths$base, cellType, trajectory)

message("Starting switch gene perturbation logic analysis...")
message("Cell Type: ", cellType)
message("Trajectory: ", trajectory)

# =============================================================================
# Load Required Data
# =============================================================================

message("\nLoading data files...")

# Load switch gene perturbation results (single KD and OE)
switchSingleKD <- loadObject(ptPaths$switchSingleTargetsKD, config, "switch single KD results")
switchSingleOE <- loadObject(ptPaths$switchSingleTargetsOE, config, "switch single OE results")

# Load GeneSwitches results for metadata
geneSwitches <- loadObject(ptPaths$geneSwitchesTsv, config, "GeneSwitches results")

message("Data loaded successfully")
message("  - Switch KD perturbations: ", nrow(switchSingleKD))
message("  - Switch OE perturbations: ", nrow(switchSingleOE))
message("  - GeneSwitches metadata: ", nrow(geneSwitches))

# =============================================================================
# PART 1: Create Unified Switch Gene Report
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PART 1: Creating Unified Switch Gene Report")
message(paste(rep("=", 60), collapse = ""))

# Get unique switch genes (intersection of KD and OE results)
kdGenes <- switchSingleKD$Gene[switchSingleKD$Success]
oeGenes <- switchSingleOE$Gene[switchSingleOE$Success]
allSwitchGenes <- intersect(kdGenes, oeGenes)

message("Building comprehensive report for ", length(allSwitchGenes), " switch genes...")

# Initialize unified report
switchReport <- data.frame(
  Gene = allSwitchGenes,
  SwitchDirection = NA_character_,
  Pseudotime = NA_real_,
  KD_Delta = NA_real_,
  OE_Delta = NA_real_,
  Category = NA_character_,
  stringsAsFactors = FALSE
)

# Fill in data for each gene
for (i in seq_along(allSwitchGenes)) {
  gene <- allSwitchGenes[i]
  
  # Get GeneSwitches metadata
  gsRow <- geneSwitches[geneSwitches$geneId == gene, ]
  if (nrow(gsRow) > 0) {
    switchReport$SwitchDirection[i] <- gsRow$direction[1]
    switchReport$Pseudotime[i] <- gsRow$pseudotime[1]
  }
  
  # Get KD results
  kdRow <- switchSingleKD[switchSingleKD$Gene == gene, ]
  if (nrow(kdRow) > 0) {
    switchReport$KD_Delta[i] <- kdRow$Delta[1]
  }
  
  # Get OE results
  oeRow <- switchSingleOE[switchSingleOE$Gene == gene, ]
  if (nrow(oeRow) > 0) {
    switchReport$OE_Delta[i] <- oeRow$Delta[1]
  }
}

message("Part 1 completed: Basic report structure created")

# =============================================================================
# PART 2: Classify Switch Genes into Categories
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PART 2: Classifying Switch Genes by Perturbation Logic")
message(paste(rep("=", 60), collapse = ""))

# Classify each gene based on perturbation outcomes
for (i in seq_len(nrow(switchReport))) {
  switchDir <- switchReport$SwitchDirection[i]
  kdDelta <- switchReport$KD_Delta[i]
  oeDelta <- switchReport$OE_Delta[i]
  
  if (is.na(switchDir) || is.na(kdDelta) || is.na(oeDelta)) {
    switchReport$Category[i] <- "Incomplete Data"
    next
  }
  
  # Determine reinforcing and opposing deltas
  if (switchDir == "up") {
    reinforcingDelta <- oeDelta  # OE reinforces upward switch
    opposingDelta <- kdDelta     # KD opposes upward switch
  } else {  # switchDir == "down"
    reinforcingDelta <- kdDelta  # KD reinforces downward switch
    opposingDelta <- oeDelta     # OE opposes downward switch
  }
  
  # Use threshold for determining "worsens" or "improves"
  threshold <- 0.001
  reinforcingWorsens <- reinforcingDelta > threshold
  opposingImproves <- opposingDelta < -threshold
  opposingWorsens <- opposingDelta > threshold
  reinforcingImproves <- reinforcingDelta < -threshold
  
  # Classify into categories
  if (reinforcingWorsens && opposingImproves) {
    # Pattern 1 & 2: Pro-aging driver (both directions confirm it)
    switchReport$Category[i] <- "Pro-aging Driver"
    
  } else if (opposingWorsens && reinforcingImproves) {
    # Pattern 3 & 4: Protective/beneficial (both directions confirm it)
    switchReport$Category[i] <- "Protective"
    
  } else if (reinforcingWorsens && !opposingImproves) {
    # Reinforcing is bad, but opposing doesn't help
    switchReport$Category[i] <- "Pro-aging (weak)"
    
  } else if (opposingImproves && !reinforcingWorsens) {
    # Opposing helps, but reinforcing doesn't hurt
    switchReport$Category[i] <- "Rejuvenation Target (weak)"
    
  } else if (opposingWorsens && !reinforcingImproves) {
    # Opposing is bad, but reinforcing doesn't help
    switchReport$Category[i] <- "Protective (weak)"
    
  } else if (reinforcingImproves && !opposingWorsens) {
    # Reinforcing helps, but opposing doesn't hurt
    switchReport$Category[i] <- "Beneficial Adaptation (weak)"
    
  } else {
    # Neither direction has clear effect
    switchReport$Category[i] <- "Ambiguous"
  }
}

# Count categories
categoryCounts <- table(switchReport$Category)
message("\nCategory Distribution:")
for (cat in names(categoryCounts)) {
  message("  ", cat, ": ", categoryCounts[cat])
}

message("\nPart 2 completed: All genes classified")

# =============================================================================
# PART 3: Sort by Pseudotime
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PART 3: Sorting by Pseudotime")
message(paste(rep("=", 60), collapse = ""))

# Sort by pseudotime (when the switch occurs)
switchReport <- switchReport[order(switchReport$Pseudotime), ]

message("Part 3 completed: Report sorted by pseudotime")

# =============================================================================
# Save Results
# =============================================================================

if (config$saveResults) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Saving Results")
  message(paste(rep("=", 60), collapse = ""))
  
  # Save unified switch gene report (TSV only)
  write.table(switchReport, ptPaths$switchPerturbationAnalysisTsv, 
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Report saved: ", basename(ptPaths$switchPerturbationAnalysisTsv))
  
  message("\nResults saved successfully")
}

# =============================================================================
# Final Summary
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SWITCH GENE PERTURBATION LOGIC ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("Cell Type: ", cellType)
message("Trajectory: ", trajectory)
message("Total switch genes analyzed: ", nrow(switchReport))
message("")
message("CATEGORY BREAKDOWN:")
for (cat in names(categoryCounts)) {
  message("  ", cat, ": ", categoryCounts[cat])
}
message("")

# Highlight top targets from key categories
proAgingDrivers <- switchReport[switchReport$Category == "Pro-aging Driver", ]
protectiveGenes <- switchReport[switchReport$Category == "Protective", ]

if (nrow(proAgingDrivers) > 0) {
  message("PRO-AGING DRIVERS (", nrow(proAgingDrivers), " genes):")
  topDrivers <- head(proAgingDrivers[, c("Gene", "SwitchDirection", "Pseudotime", "KD_Delta", "OE_Delta")], 5)
  for (i in seq_len(nrow(topDrivers))) {
    message("  ", topDrivers$Gene[i], 
            " (", topDrivers$SwitchDirection[i], " @ t=", round(topDrivers$Pseudotime[i], 2), 
            ") - KD Δ=", round(topDrivers$KD_Delta[i], 4),
            ", OE Δ=", round(topDrivers$OE_Delta[i], 4))
  }
  message("")
}

if (nrow(protectiveGenes) > 0) {
  message("PROTECTIVE GENES (", nrow(protectiveGenes), " genes):")
  topProtective <- head(protectiveGenes[, c("Gene", "SwitchDirection", "Pseudotime", "KD_Delta", "OE_Delta")], 5)
  for (i in seq_len(nrow(topProtective))) {
    message("  ", topProtective$Gene[i],
            " (", topProtective$SwitchDirection[i], " @ t=", round(topProtective$Pseudotime[i], 2),
            ") - KD Δ=", round(topProtective$KD_Delta[i], 4),
            ", OE Δ=", round(topProtective$OE_Delta[i], 4))
  }
  message("")
}
message(paste(rep("=", 60), collapse = ""))
message("Done!")
