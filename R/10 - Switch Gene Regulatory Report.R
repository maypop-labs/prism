# =============================================================================
# 10 - Switch Gene Regulatory Report
# Purpose: Extract switch gene regulatory relationships from GeneSwitches and
#          Boolean rules to create a comprehensive regulatory network report
# =============================================================================

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
ensureProjectDirectories(paths)
clearConsole()

# --- Load Data ---
message("Loading GeneSwitches results...")
switchGenes <- loadObject(ptPaths$geneSwitches, config, "switch genes")
switchDf <- as.data.frame(switchGenes, stringsAsFactors = FALSE)

message("Loading Boolean rules...")
boolRules <- loadObject(ptPaths$booleanRules, config, "Boolean rules")

message("Loading gene name mapping...")
geneMap <- loadObject(ptPaths$geneMap, config, "gene mapping")

# --- Create reverse lookup for sanitized names ---
reverseLookup <- setNames(geneMap$OriginalName, geneMap$SanitizedName)

# --- Create set of switch gene names for filtering ---
# Include both original and sanitized names
switchGeneNamesOriginal <- unique(switchDf$geneID)
switchGeneNamesSanitized <- unique(geneMap$SanitizedName[geneMap$OriginalName %in% switchDf$geneID])
allSwitchGeneNames <- unique(c(switchGeneNamesOriginal, switchGeneNamesSanitized))

# --- Parse Boolean rules to extract regulatory relationships ---
message("Parsing Boolean rules for regulatory relationships...")

# Initialize results list
regulatoryInfo <- list()

for (targetGene in names(boolRules)) {
  rule <- boolRules[[targetGene]]
  
  # Extract rule string
  ruleString <- rule$rule %||% ""
  
  if (nchar(ruleString) == 0 || ruleString == targetGene) {
    # Self-activation or empty rule - no regulators
    regulatoryInfo[[targetGene]] <- list(
      regulators = character(0),
      regTypes = character(0)
    )
    next
  }
  
  # Extract all gene names from the rule (both with and without !)
  # Pattern matches valid R identifiers that might be genes
  allMatches <- gregexpr("[A-Za-z][A-Za-z0-9_.]*", ruleString)
  allGenes <- regmatches(ruleString, allMatches)[[1]]
  
  # Filter out the target gene itself and logical operators
  logicalOps <- c("TRUE", "FALSE", "NA", "NULL")
  potentialRegulators <- setdiff(allGenes, c(targetGene, logicalOps))
  
  if (length(potentialRegulators) == 0) {
    regulatoryInfo[[targetGene]] <- list(
      regulators = character(0),
      regTypes = character(0)
    )
    next
  }
  
  # Determine activation vs repression for each regulator
  regulators <- character()
  regTypes <- character()
  
  for (reg in potentialRegulators) {
    # Check if regulator appears with negation
    # Look for patterns like !reg or !(...)reg(...)
    negPattern <- paste0("!", reg, "(?![A-Za-z0-9_])")
    hasNegation <- grepl(negPattern, ruleString, perl = TRUE)
    
    regulators <- c(regulators, reg)
    regTypes <- c(regTypes, ifelse(hasNegation, "Repression", "Activation"))
  }
  
  regulatoryInfo[[targetGene]] <- list(
    regulators = regulators,
    regTypes = regTypes
  )
}

# --- Build comprehensive report ---
message("Building switch gene regulatory report...")

# Create a data frame with one row per switch gene
reportData <- data.frame(
  switchGene = character(),
  switchGene_original = character(),
  direction = character(),
  pseudotime = numeric(),
  fdr = numeric(),
  pseudoR2 = numeric(),
  regulatesGenes = character(),
  regulationTypes = character(),
  n_targets = integer(),
  stringsAsFactors = FALSE
)

# Process each switch gene
for (i in seq_len(nrow(switchDf))) {
  geneId <- switchDf$geneID[i]
  
  # Try to find the sanitized version
  sanitizedGene <- NA
  if (geneId %in% geneMap$OriginalName) {
    sanitizedGene <- geneMap$SanitizedName[geneMap$OriginalName == geneId]
  } else {
    sanitizedGene <- geneId  # Already sanitized or not in map
  }
  
  # Check if this switch gene acts as a regulator in any rules
  regulatedGenes <- character()
  regulationTypes <- character()
  
  for (targetGene in names(regulatoryInfo)) {
    regInfo <- regulatoryInfo[[targetGene]]
    
    # Check if our switch gene is a regulator of this target
    if (sanitizedGene %in% regInfo$regulators) {
      # Convert target back to original name if possible
      targetOriginal <- reverseLookup[[targetGene]] %||% targetGene
      
      # FILTER: Only include if target is also a switch gene
      if (targetOriginal %in% allSwitchGeneNames || targetGene %in% allSwitchGeneNames) {
        # Get regulation type
        idx <- which(regInfo$regulators == sanitizedGene)
        regType <- regInfo$regTypes[idx]
        
        regulatedGenes <- c(regulatedGenes, targetOriginal)
        regulationTypes <- c(regulationTypes, regType)
      }
    }
  }
  
  # Combine targets and types into single strings
  targetsString <- if (length(regulatedGenes) == 0) "" else paste(regulatedGenes, collapse = "; ")
  typesString <- if (length(regulationTypes) == 0) "" else paste(regulationTypes, collapse = "; ")
  
  # Add to report
  reportData <- rbind(reportData, data.frame(
    switchGene = sanitizedGene,
    switchGene_original = geneId,
    direction = switchDf$direction[i],
    pseudotime = switchDf$switch_at_time[i],
    fdr = switchDf$FDR[i],
    pseudoR2 = switchDf$pseudoR2s[i],
    regulatesGenes = targetsString,
    regulationTypes = typesString,
    n_targets = length(regulatedGenes),
    stringsAsFactors = FALSE
  ))
}

# Sort by pseudotime (ascending order)
reportData <- reportData[order(reportData$pseudotime), ]
rownames(reportData) <- NULL

# --- Generate summary statistics ---
message("\nSwitch Gene Regulatory Network Summary:")
message("  Total switch genes: ", nrow(reportData))
message("  Switch genes that regulate other switch genes: ", sum(reportData$n_targets > 0))
message("  Switch genes with no switch gene targets: ", sum(reportData$n_targets == 0))
message("  Total switch-to-switch regulatory relationships: ", sum(reportData$n_targets))

if (sum(reportData$n_targets > 0) > 0) {
  message("  Mean switch gene targets per regulator: ", 
          round(mean(reportData$n_targets[reportData$n_targets > 0]), 2))
}

# Count regulation types
allTypes <- unlist(strsplit(reportData$regulationTypes[reportData$regulationTypes != ""], "; "))
if (length(allTypes) > 0) {
  typeTable <- table(allTypes)
  message("\nRegulation types (switch-to-switch only):")
  for (type in names(typeTable)) {
    message("  ", type, ": ", typeTable[type])
  }
}

# --- Save report ---
if (config$saveResults) {
  outputPath <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_switch_regulatory_network.tsv"))
  write.table(reportData, outputPath, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\nReport saved to: ", outputPath)
  
  # Also save a simplified version with just key columns
  simplifiedReport <- reportData[, c("switchGene_original", "direction", "pseudotime", 
                                      "regulatesGenes", "regulationTypes", "n_targets")]
  simplifiedPath <- file.path(paths$base$tsv, paste0(cellType, "_", trajectory, "_switch_regulatory_network_simple.tsv"))
  write.table(simplifiedReport, simplifiedPath, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Simplified report saved to: ", simplifiedPath)
}

# =============================================================================
# Network Visualization
# =============================================================================

message("\n--- Creating Network Visualization ---")

# Only create visualization if there are regulatory relationships
if (sum(reportData$n_targets > 0) > 0) {
  
  # --- Build edge list from reportData ---
  edgeList <- data.frame(
    from = character(),
    to = character(),
    regType = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(reportData))) {
    if (reportData$regulatesGenes[i] != "") {
      # Parse the semicolon-separated targets and types
      targets <- strsplit(reportData$regulatesGenes[i], "; ")[[1]]
      types <- strsplit(reportData$regulationTypes[i], "; ")[[1]]
      
      for (j in seq_along(targets)) {
        edgeList <- rbind(edgeList, data.frame(
          from = reportData$switchGene_original[i],
          to = targets[j],
          regType = types[j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # --- Create node attributes ---
  nodeAttrs <- reportData[, c("switchGene_original", "direction", "pseudotime", "n_targets")]
  colnames(nodeAttrs) <- c("name", "direction", "pseudotime", "n_targets")
  
  # --- Build igraph network ---
  message("Building network graph...")
  net <- graph_from_data_frame(d = edgeList, vertices = nodeAttrs, directed = TRUE)
  
  # --- Create custom layout with pseudotime on Y-axis ---
  # Start with a constrained force-directed layout
  set.seed(config$randomSeed)  # Use config seed for reproducibility
  
  # Get number of nodes
  nNodes <- vcount(net)
  
  # Create initial layout with Y-axis constrained to pseudotime
  # X positions will be optimized to minimize edge crossings
  nodeData <- as.data.frame(vertex_attr(net))
  
  # Strategy: Use stress layout but constrain Y to pseudotime
  # First, get an unconstrained layout
  layoutInit <- layout_with_fr(net, niter = 500)
  
  # Replace Y coordinates with pseudotime (scaled and inverted)
  # Invert so early pseudotime (young) is at top
  maxPseudo <- max(nodeData$pseudotime)
  minPseudo <- min(nodeData$pseudotime)
  yCoords <- -(nodeData$pseudotime - minPseudo) / (maxPseudo - minPseudo)
  
  # For X coordinates, use the layout's X but scale to reasonable range
  xCoords <- layoutInit[, 1]
  # Scale X to be narrower (creates vertical column effect)
  xCoords <- scale(xCoords) * 0.3
  
  # Combine into final layout
  layoutFinal <- cbind(xCoords, yCoords)
  
  # --- Create visualization using ggraph ---
  message("Creating visualization...")
  
  # Set up colors
  directionColors <- c("up" = "#E74C3C", "down" = "#3498DB")  # Red for up, blue for down
  regTypeColors <- c("Activation" = "#27AE60", "Repression" = "#E67E22")  # Green for activation, orange for repression
  
  # Calculate node size scaling based on config pointSize
  baseNodeSize <- config$pointSize * 2  # Base multiplier for node size
  nodeSizeRange <- c(baseNodeSize, baseNodeSize * 5)
  
  # Create the plot
  p <- ggraph(net, layout = layoutFinal) +
    # Draw edges first (so they're behind nodes)
    geom_edge_link(
      aes(color = regType),
      arrow = arrow(length = unit(3, 'mm'), type = "closed"),
      end_cap = circle(3, 'mm'),
      start_cap = circle(3, 'mm'),
      alpha = config$plotAlpha * 0.8,  # Slightly more transparent for edges
      width = 0.8
    ) +
    # Draw nodes
    geom_node_point(
      aes(size = n_targets + 1, fill = direction),
      shape = 21,
      color = "black",
      stroke = 0.5,
      alpha = config$plotAlpha
    ) +
    # Add gene labels
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 2.5,
      fontface = "italic",
      max.overlaps = 20
    ) +
    # Color scales
    scale_edge_color_manual(
      values = regTypeColors,
      name = "Regulation Type"
    ) +
    scale_fill_manual(
      values = directionColors,
      name = "Switch Direction",
      labels = c("down" = "Downregulated", "up" = "Upregulated")
    ) +
    scale_size_continuous(
      name = "Target Count",
      range = nodeSizeRange,
      breaks = c(1, 3, 5, 10)
    ) +
    # Labels and theme
    labs(
      title = paste0("Switch Gene Regulatory Network: ", cellType, " - ", trajectory),
      subtitle = "Temporal cascade of switch genes (top = early aging, bottom = late aging)",
      caption = paste0(
        "Nodes: switch genes (size = # of switch gene targets)\n",
        "Edges: regulatory relationships between switch genes\n",
        "Y-axis position: pseudotime of switch occurrence"
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray50"),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  # --- Save plot ---
  if (config$saveResults) {
    outputPath <- file.path(paths$base$plots, 
                           paste0(cellType, "_", trajectory, "_switch_regulatory_network.png"))
    
    # Determine appropriate dimensions based on network size
    # Use config default width, but scale height with network size
    plotWidth <- config$figWidth
    plotHeight <- max(config$figHeight, nNodes * 0.4)  # Scale with nodes
    plotHeight <- min(plotHeight, 20)   # Cap at 20 inches
    
    ggsave(
      filename = outputPath,
      plot = p,
      width = plotWidth,
      height = plotHeight,
      dpi = config$figDPI,
      bg = "white"
    )
    
    message("Network visualization saved to: ", outputPath)
  }
  
} else {
  message("No regulatory relationships found between switch genes - skipping visualization")
}

message("\nDone!")
