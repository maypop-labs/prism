# =============================================================================
# booleanManager.R
# Purpose: Boolean network construction, rule inference, and optimization
# Dependencies: BoolNet, foreach, doParallel
# =============================================================================

# =============================================================================
# Synthesize Boolean Rules from SCENIC with Sign Annotations
# =============================================================================

#' Synthesize Boolean rule for a target gene based on signed SCENIC edges
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector (1 for activation, -1 for repression)
#' @return Character string representing a simple Boolean rule
synthesizeBooleanRule <- function(targetGene, regulators, regulatorSigns) {
  if (length(regulators) == 0) {
    return(paste0(targetGene, ", ", targetGene))
  }
  
  terms <- mapply(function(reg, sign) {
    if (sign > 0) {
      return(reg)
    } else if (sign < 0) {
      return(paste0("!", reg))
    } else {
      return(NULL)
    }
  }, regulators, regulatorSigns[regulators], SIMPLIFY = TRUE)
  
  # Combine with AND by default (could be OR if specified)
  ruleLogic <- paste(terms, collapse = " & ")
  return(paste0(targetGene, ", ", ruleLogic))
}

# =============================================================================
# Batch Synthesis from SCENIC Edge List
# =============================================================================

#' Synthesize Boolean rules for all target genes using SCENIC edges
#'
#' @param edges Data frame with columns TF, Target, regType
#' @param matBin Binary expression matrix
#' @param maxRegulators Maximum number of regulators per target
#' @return Named list of rule strings
synthesizeBooleanRulesBatch <- function(edges, matBin, maxRegulators = 5) {
  rules <- list()
  targets <- unique(edges$Target)
  total <- length(targets)
  
  message("Synthesizing rules for ", total, " target genes...")
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in seq_along(targets)) {
    gene <- targets[i]
    setTxtProgressBar(pb, i)
    
    subEdges <- edges[edges$Target == gene, ]
    if (nrow(subEdges) == 0) next
    
    topEdges <- head(subEdges[order(-abs(subEdges$corr)), ], maxRegulators)
    regulators <- topEdges$TF
    signs <- setNames(ifelse(topEdges$regType == "activation", 1, -1), topEdges$TF)
    
    ruleStr <- synthesizeBooleanRule(gene, regulators, signs)
    score <- scoreBooleanRule(ruleStr, matBin)
    
    rules[[gene]] <- list(
      rule = ruleStr,
      regulators = regulators,
      score = score,
      nRegulators = length(regulators),
      methodUsed = "template",
      biologicallyPlausible = NA
    )
  }
  
  close(pb)
  return(rules)
}

# =============================================================================
# Enhanced Boolean Rule Synthesis with Multiple Templates
# =============================================================================

#' Generate multiple Boolean rule templates for a target gene
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector (1 for activation, -1 for repression)
#' @return List of rule templates with descriptions
generateBooleanTemplates <- function(targetGene, regulators, regulatorSigns) {
  
  if (length(regulators) == 0) {
    return(list(
      self_activation = list(
        rule = paste0(targetGene, ", ", targetGene),
        description = "Self-activation (default state)"
      )
    ))
  }
  
  # Prepare terms with signs
  activators <- regulators[regulatorSigns[regulators] > 0]
  repressors <- regulators[regulatorSigns[regulators] < 0]
  
  templates <- list()
  
  # 1. SINGLE REGULATOR TEMPLATES (if only 1 regulator)
  if (length(regulators) == 1) {
    reg <- regulators[1]
    if (regulatorSigns[reg] > 0) {
      templates$single_activator <- list(
        rule = paste0(targetGene, ", ", reg),
        description = paste("Single activator:", reg)
      )
    } else {
      templates$single_repressor <- list(
        rule = paste0(targetGene, ", !", reg),
        description = paste("Single repressor:", reg)
      )
    }
    return(templates)
  }
  
  # 2. COOPERATIVE ACTIVATION (AND logic)
  if (length(activators) >= 2) {
    andTerms <- c(activators, if(length(repressors) > 0) paste0("!", repressors) else character(0))
    templates$cooperative_activation <- list(
      rule = paste0(targetGene, ", ", paste(andTerms, collapse = " & ")),
      description = "Cooperative activation (all activators + no repressors)"
    )
  }
  
  # 3. REDUNDANT ACTIVATION (OR logic for activators)
  if (length(activators) >= 2) {
    # Any activator works, but repressors can block
    if (length(repressors) > 0) {
      orActivators <- paste0("(", paste(activators, collapse = " | "), ")")
      andRepressors <- paste0("!(", paste(repressors, collapse = " | "), ")")
      ruleLogic <- paste(orActivators, andRepressors, sep = " & ")
    } else {
      ruleLogic <- paste(activators, collapse = " | ")
    }
    
    templates$redundant_activation <- list(
      rule = paste0(targetGene, ", ", ruleLogic),
      description = "Redundant activation (any activator works)"
    )
  }
  
  # 4. BALANCED ACTIVATION (mix of AND and OR)
  if (length(activators) >= 2 && length(repressors) >= 1) {
    # Need at least one activator AND no repressors
    orActivators <- paste0("(", paste(activators, collapse = " | "), ")")
    andRepressors <- if(length(repressors) > 0) paste0("!", paste(repressors, collapse = " & !")) else ""
    
    templates$balanced_regulation <- list(
      rule = paste0(targetGene, ", ", orActivators, " & ", andRepressors),
      description = "Balanced: any activator + all repressors off"
    )
  }
  
  if (length(activators) >= 2 && length(repressors) >= 1) {
    orActivators <- paste0("(", paste(activators, collapse = " | "), ")")
    andRepressors <- paste0("!", paste(repressors, collapse = " & !"))
    
    templates$balanced_regulation <- list(
      rule = paste0(targetGene, ", ", orActivators, " & ", andRepressors),
      description = "Balanced: any activator + all repressors off"
    )
  }
  
  # 5. MAJORITY RULE (for 3+ regulators)
  if (length(regulators) >= 3) {
    # Simple majority: more activators than repressors
    majorityTerms <- character()
    
    # All possible combinations where activators outnumber active repressors
    if (length(activators) >= 2) {
      # At least 2 activators on
      activatorPairs <- combn(activators, 2, simplify = FALSE)
      for (pair in activatorPairs) {
        pairTerm <- paste(pair, collapse = " & ")
        majorityTerms <- c(majorityTerms, pairTerm)
      }
      
      # If we have 3+ activators, add triple combinations
      if (length(activators) >= 3) {
        allActivators <- paste(activators, collapse = " & ")
        majorityTerms <- c(majorityTerms, allActivators)
      }
    }
    
    if (length(majorityTerms) > 0) {
      templates$majority_rule <- list(
        rule = paste0(targetGene, ", ", paste(unique(majorityTerms), collapse = " | ")),
        description = "Majority rule (multiple activator combinations)"
      )
    }
  }
  
  # 6. SIMPLE AND (default fallback)
  allTerms <- c(activators, if(length(repressors) > 0) paste0("!", repressors) else character(0))
  templates$simple_and <- list(
    rule = paste0(targetGene, ", ", paste(allTerms, collapse = " & ")),
    description = "Simple AND (all activators + no repressors)"
  )
  
  # 7. DOMINANT ACTIVATOR (strongest activator alone)
  if (length(activators) >= 1) {
    # Assume first activator is strongest (could sort by SCENIC score)
    dominantActivator <- activators[1]
    templates$dominant_activator <- list(
      rule = paste0(targetGene, ", ", dominantActivator),
      description = paste("Dominant activator:", dominantActivator)
    )
  }
  
  return(templates)
}

#' Test multiple templates and select the best one
#'
#' @param targetGene Character name of the target gene
#' @param regulators Character vector of regulator gene names
#' @param regulatorSigns Named numeric vector of regulatory signs
#' @param matBin Binary expression matrix
#' @return List with best rule, score, and all tested templates
synthesizeBestBooleanRule <- function(targetGene, regulators, regulatorSigns, matBin) {
  
  # Generate all possible templates
  templates <- generateBooleanTemplates(targetGene, regulators, regulatorSigns)
  
  # Score each template
  templateScores <- list()
  bestScore <- -1
  bestTemplate <- NULL
  bestRuleName <- NULL
  
  for (templateName in names(templates)) {
    template <- templates[[templateName]]
    ruleStr <- template$rule
    
    tryCatch({
      score <- scoreBooleanRule(ruleStr, matBin)
      templateScores[[templateName]] <- list(
        rule = ruleStr,
        score = score,
        description = template$description
      )
      
      if (score > bestScore) {
        bestScore <- score
        bestTemplate <- template
        bestRuleName <- templateName
      }
    }, error = function(e) {
      # Skip templates that fail to evaluate
      templateScores[[templateName]] <- list(
        rule = ruleStr,
        score = 0,
        description = paste(template$description, "(FAILED)")
      )
    })
  }
  
  return(list(
    bestRule = bestTemplate$rule,
    bestScore = bestScore,
    bestTemplateName = bestRuleName,
    bestDescription = bestTemplate$description,
    allTemplates = templateScores,
    regulators = regulators,
    nRegulators = length(regulators)
  ))
}

#' Enhanced batch synthesis with template variety
#'
#' @param edges Data frame with columns TF, Target, regType
#' @param matBin Binary expression matrix
#' @param maxRegulators Maximum number of regulators per target
#' @return Named list of rule objects with template information
synthesizeBooleanRulesBatchEnhanced <- function(edges, matBin, maxRegulators = 5) {
  
  # DEBUGGING: Check inputs
  message("=== DEBUG START ===")
  message("Edges dimensions: ", nrow(edges), " x ", ncol(edges))
  message("MatBin dimensions: ", nrow(matBin), " x ", ncol(matBin))
  message("Max regulators: ", maxRegulators)
  
  targets <- unique(edges$Target)
  total <- length(targets)
  message("Unique targets found: ", total)
  
  if (total == 0) {
    message("ERROR: No targets found!")
    return(list())
  }
  
  # Check first few targets
  message("First 5 targets: ", paste(head(targets, 5), collapse = ", "))
  
  rules <- list()
  message("Starting synthesis loop...")
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in seq_along(targets)) {
    gene <- targets[i]
    setTxtProgressBar(pb, i)
    
    # Get edges for this gene
    subEdges <- edges[edges$Target == gene, ]

    if (nrow(subEdges) == 0) {
      message("No edges found for ", gene, " - skipping")
      next
    }
    
    # Get regulators
    topEdges <- head(subEdges[order(-abs(subEdges$corr)), ], maxRegulators)
    regulators <- topEdges$TF
    
    # Check if gene and regulators are in matrix
    geneInMatrix <- gene %in% rownames(matBin)
    regsInMatrix <- regulators %in% rownames(matBin)
    
    if (!geneInMatrix) {
      message("Gene not in matrix - skipping")
      next
    }
    
    if (!any(regsInMatrix)) {
      message("No regulators in matrix - skipping")
      next
    }
    
    # Create signs
    signs <- setNames(ifelse(topEdges$regType == "Activation", 1, -1), topEdges$TF)
    
    # Try to synthesize rule
    tryCatch({
      ruleResult <- synthesizeBestBooleanRule(gene, regulators, signs, matBin)
      
      rules[[gene]] <- list(
        rule = ruleResult$bestRule,
        regulators = regulators,
        bestScore = ruleResult$bestScore,
        score = ruleResult$bestScore,
        nRegulators = length(regulators),
        templateUsed = ruleResult$bestTemplateName,
        templateDescription = ruleResult$bestDescription,
        methodUsed = "template_variety",
        biologicallyPlausible = TRUE
      )
      
    }, error = function(e) {
      message("ERROR in synthesizeBestBooleanRule for ", gene, ": ", e$message)
    })
    
  }
  
  close(pb)
  return(rules)
}

#' Helper function to score Boolean rules (same as before but with better error handling)
#'
#' @param ruleStr BoolNet-style rule string
#' @param matBin Binary expression matrix
#' @return Accuracy score (0 to 1)
scoreBooleanRule <- function(ruleStr, matBin) {
  
  ruleParts <- strsplit(ruleStr, ",\\s*")[[1]]
  if (length(ruleParts) != 2) {
    warning("Invalid rule format: ", ruleStr)
    return(0)
  }
  
  targetGene <- ruleParts[1]
  logic <- ruleParts[2]
  
  if (!targetGene %in% rownames(matBin)) {
    warning("Target gene not found in matrix: ", targetGene)
    return(0)
  }
  
  n <- ncol(matBin)
  predictions <- logical(n)
  actual <- as.logical(matBin[targetGene, ])
  
  # Check if all required genes are in the matrix
  requiredGenes <- unique(gsub("[!&|() ]", "", unlist(strsplit(logic, "\\b"))))
  requiredGenes <- requiredGenes[requiredGenes != "" & requiredGenes %in% rownames(matBin)]
  
  if (length(requiredGenes) == 0) {
    warning("No valid regulator genes found for rule: ", ruleStr)
    return(0)
  }
  
  for (i in 1:n) {
    env <- as.list(as.logical(matBin[, i]))
    names(env) <- rownames(matBin)
    
    predictions[i] <- tryCatch({
      result <- eval(parse(text = logic), envir = env)
      if (is.logical(result) && length(result) == 1 && !is.na(result)) {
        result
      } else {
        FALSE
      }
    }, error = function(e) {
      FALSE
    })
  }
  
  return(mean(predictions == actual, na.rm = TRUE))
}

# =============================================================================
# Helper: Convert Synthesized Rules to BoolNet Format File
# =============================================================================

#' Write BoolNet-compatible rules to file
#'
#' @param ruleList List of synthesized rules
#' @param outputPath Path to write the network file
#' @export
writeBoolNetFile <- function(ruleList, outputPath) {
  lines <- sapply(ruleList, function(x) x$rule)
  writeLines(lines, con = outputPath)
}

# =============================================================================
# Boolean Rule Analysis and Visualization Functions
# =============================================================================

#' Generate comprehensive Boolean rule analysis report
#'
#' Creates visualizations using config parameters and proper white backgrounds
#'
#' @param boolRules List of Boolean rules from main script
#' @param edges Edge list data frame with regulatory relationships
#' @param paths Paths object with directory structure
#' @param cellType Cell type name for file naming
#' @param trajectory Trajectory name for file naming
#' @param config Configuration object with plot parameters
#' @export
generateBooleanRuleReport <- function(boolRules, edges, paths, cellType, trajectory, config) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(ComplexHeatmap)
  library(igraph)
  
  if (!dir.exists(paths$base$plots)) dir.create(paths$base$plots, recursive = TRUE)
  if (!dir.exists(paths$base$txt)) dir.create(paths$base$txt, recursive = TRUE)
  
  # Extract summary statistics
  ruleStats <- extractRuleStatistics(boolRules)
  
  # Standard plot theme using config parameters
  standardTheme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  # 1. Rule Quality Distribution
  p1 <- plotRuleQualityDistribution(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_rule_quality_distribution.png"), 
         p1, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 2. Method Usage Summary
  p2 <- plotMethodUsageSummary(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_method_usage_summary.png"), 
         p2, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 3. Biological Plausibility Analysis
  p3 <- plotBiologicalPlausibility(ruleStats, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_biological_plausibility.png"), 
         p3, width = config$figWidth, height = config$figHeight, dpi = config$figDPI, bg = "white")
  
  # 4. Network Complexity Heatmap (FIXED VERSION)
  if (length(boolRules) > 0) {
    p4 <- plotNetworkComplexityHeatmapFixed(boolRules, edges)
    png(paste0(paths$base$plots, cellType, "_", trajectory, "_network_complexity_heatmap.png"), 
        width = config$figWidth * config$figDPI, height = config$figHeight * config$figDPI, 
        res = config$figDPI, bg = "white")
    draw(p4)
    dev.off()
  }
  
  # 5. Rule Logic Summary Network
  p5 <- plotRuleLogicNetwork(boolRules, edges, standardTheme)
  ggsave(paste0(paths$base$plots, cellType, "_", trajectory, "_rule_logic_network.png"), 
         p5, width = config$figWidth * 1.5, height = config$figHeight * 1.2, dpi = config$figDPI, bg = "white")
  
  # 6. Generate text report
  generateTextReport(ruleStats, boolRules, paste0(paths$base$txt, cellType, "_", trajectory, "_boolean_rules_report.txt"))
  
  message("Boolean rule analysis complete. Results saved to: ", paths$base$plots)
}

#' Plot distribution of rule quality scores
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotRuleQualityDistribution <- function(ruleStats, standardTheme) {
  
  ggplot(ruleStats, aes(x = bestScore)) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.95, color = "darkgreen", linetype = "dashed", size = 1) +
    labs(
      title = "Boolean Rule Quality Distribution",
      subtitle = paste0("n = ", nrow(ruleStats), " genes with inferred rules"),
      x = "Rule Quality Score (Accuracy)",
      y = "Number of Genes",
      caption = "Red line: minimum threshold (0.8), Green line: excellent threshold (0.95)"
    ) +
    standardTheme +
    annotate("text", x = 0.85, y = Inf, label = "Good Rules", vjust = 1.5, color = "darkgreen", size = 3) +
    annotate("text", x = 0.4, y = Inf, label = "Poor Rules", vjust = 1.5, color = "red", size = 3)
}

#' Plot method usage summary
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotMethodUsageSummary <- function(ruleStats, standardTheme) {
  
  methodCounts <- ruleStats %>%
    dplyr::group_by(methodUsed) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::mutate(
      percentage = round(100 * count / sum(count), 1),
      label = paste0(percentage, "%")
    )
  
  ggplot(methodCounts, aes(x = "", y = count, fill = methodUsed)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(
      title = "Boolean Rule Inference Method Usage",
      subtitle = "Which methods were used to find the best rules?",
      fill = "Method"
    ) +
    standardTheme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
}

#' Plot biological plausibility analysis
#' @param ruleStats Data frame with rule statistics
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotBiologicalPlausibility <- function(ruleStats, standardTheme) {
  
  plausibilityData <- ruleStats %>%
    filter(!is.na(biologicallyPlausible)) %>%
    mutate(
      plausibilityLabel = ifelse(biologicallyPlausible, "Biologically Plausible", "Data-Driven Only")
    )
  
  if (nrow(plausibilityData) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No biological plausibility data available", size = 5) +
             labs(title = "Biological Plausibility Analysis", subtitle = "No data available") +
             standardTheme)
  }
  
  ggplot(plausibilityData, aes(x = plausibilityLabel, y = bestScore, fill = plausibilityLabel)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    labs(
      title = "Rule Quality by Biological Plausibility",
      subtitle = "Do biologically plausible rules perform better?",
      x = "Rule Type",
      y = "Quality Score",
      fill = "Rule Type"
    ) +
    standardTheme +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("steelblue", "orange")) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, color = "red", fill = "white")
}

#' Create network complexity heatmap
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @return ComplexHeatmap object
plotNetworkComplexityHeatmapFixed <- function(boolRules, edges) {
  
  # FIXED: Only use genes that actually have Boolean rules and their regulators
  targetGenes <- names(boolRules)
  regulatorGenes <- unique(unlist(lapply(boolRules, function(x) x$regulators)))
  relevantGenes <- unique(c(targetGenes, regulatorGenes))
  
  # Create appropriately sized matrix
  regulatorMatrix <- matrix(0, nrow = length(regulatorGenes), ncol = length(targetGenes))
  rownames(regulatorMatrix) <- regulatorGenes
  colnames(regulatorMatrix) <- targetGenes
  
  # Fill in regulatory relationships with quality scores
  for (gene in targetGenes) {
    if (gene %in% names(boolRules)) {
      regulators <- boolRules[[gene]]$regulators
      
      # Ensure ruleQuality is a single numeric value
      ruleQuality <- as.numeric(boolRules[[gene]]$bestScore %||% boolRules[[gene]]$score %||% 0)[1]
      if (is.na(ruleQuality)) ruleQuality <- 0
      
      for (reg in regulators) {
        if (reg %in% rownames(regulatorMatrix)) {
          regulatorMatrix[reg, gene] <- ruleQuality
        }
      }
    }
  }
  
  # Remove empty rows and columns
  nonEmptyRows <- rowSums(regulatorMatrix) > 0
  nonEmptyCols <- colSums(regulatorMatrix) > 0
  
  if (sum(nonEmptyRows) > 0 && sum(nonEmptyCols) > 0) {
    regulatorMatrix <- regulatorMatrix[nonEmptyRows, nonEmptyCols, drop = FALSE]
    
    # Create properly sized heatmap
    ComplexHeatmap::Heatmap(
      regulatorMatrix,
      name = "Rule Quality",
      col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
      cluster_rows = ifelse(nrow(regulatorMatrix) > 1, TRUE, FALSE),
      cluster_columns = ifelse(ncol(regulatorMatrix) > 1, TRUE, FALSE),
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10),
      column_title = "Target Genes",
      row_title = "Regulator Genes",
      heatmap_legend_param = list(title = "Boolean Rule\nQuality Score"),
      width = unit(min(8, ncol(regulatorMatrix) * 0.8), "cm"),
      height = unit(min(8, nrow(regulatorMatrix) * 0.8), "cm")
    )
  } else {
    # Create empty heatmap with message
    emptyMatrix <- matrix(0, nrow = 1, ncol = 1)
    rownames(emptyMatrix) <- "No Data"
    colnames(emptyMatrix) <- "No Data"
    
    ComplexHeatmap::Heatmap(
      emptyMatrix,
      name = "Rule Quality",
      col = "white",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = "No regulatory relationships to display",
      heatmap_legend_param = list(title = "Boolean Rule\nQuality Score")
    )
  }
}

#' Plot rule logic network
#' @param boolRules List of Boolean rules
#' @param edges Edge list data frame
#' @param standardTheme ggplot theme to apply
#' @return ggplot object
plotRuleLogicNetwork <- function(boolRules, edges, standardTheme) {
  
  # Filter edges to only include those with Boolean rules
  networkEdges <- edges %>%
    dplyr::filter(Target %in% names(boolRules)) %>%
    dplyr::mutate(
      ruleQuality = sapply(Target, function(x) {
        if (x %in% names(boolRules)) boolRules[[x]]$bestScore else 0
      })
    )
  
  if (nrow(networkEdges) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No network data to display", size = 5) +
             labs(title = "Gene Regulatory Network", subtitle = "No edges available") +
             standardTheme)
  }
  
  # Create igraph object
  g <- igraph::graph_from_data_frame(networkEdges, directed = TRUE)
  
  # Calculate layout
  if (igraph::vcount(g) > 1) {
    layout <- igraph::layout_with_fr(g)
  } else {
    layout <- matrix(c(0, 0), nrow = 1)
  }
  
  # Convert to data frame for ggplot
  nodeList <- data.frame(
    name = igraph::V(g)$name,
    x = layout[,1],
    y = layout[,2],
    isTarget = igraph::V(g)$name %in% names(boolRules)
  )
  
  # Create edge coordinates
  edgeList <- as.data.frame(igraph::get.edgelist(g))
  names(edgeList) <- c("from", "to")
  edgeList$quality <- igraph::E(g)$ruleQuality
  
  edgeCoords <- edgeList %>%
    dplyr::left_join(nodeList, by = c("from" = "name")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::select(-isTarget) %>%
    dplyr::left_join(nodeList, by = c("to" = "name")) %>%
    dplyr::rename(x2 = x, y2 = y, isTarget = isTarget)
  
  # Plot with better styling
  ggplot() +
    geom_segment(
      data = edgeCoords,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = quality),
      arrow = arrow(length = unit(0.15, "inches")),
      alpha = 0.7,
      size = 1
    ) +
    geom_point(
      data = nodeList,
      aes(x = x, y = y, shape = isTarget, fill = isTarget),
      size = 4,
      alpha = 0.8,
      color = "black",
      stroke = 0.5
    ) +
    scale_color_gradient(
      low = "lightgray",
      high = "red",
      name = "Rule\nQuality"
    ) +
    scale_shape_manual(
      values = c("FALSE" = 21, "TRUE" = 22),
      name = "Node Type",
      labels = c("Regulator", "Target")
    ) +
    scale_fill_manual(
      values = c("FALSE" = "lightblue", "TRUE" = "orange"),
      name = "Node Type",
      labels = c("Regulator", "Target")
    ) +
    labs(
      title = "Gene Regulatory Network with Boolean Rules",
      subtitle = paste0("Network includes ", nrow(nodeList), " genes and ", nrow(edgeCoords), " regulatory relationships")
    ) +
    standardTheme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

#' Extract summary statistics from Boolean rule list
#'
#' @param boolRules List of Boolean rule objects
#' @return Data frame with per-gene statistics
extractRuleStatistics <- function(boolRules) {
  data.frame(
    gene = names(boolRules),
    nRegulators = sapply(boolRules, function(x) x$nRegulators %||% NA),
    bestScore = sapply(boolRules, function(x) x$score %||% NA),
    methodUsed = sapply(boolRules, function(x) x$methodUsed %||% "template"),
    biologicallyPlausible = sapply(boolRules, function(x) x$biologicallyPlausible %||% NA)
  )
}

#' Generate text summary report
#' @param ruleStats Data frame with rule statistics
#' @param boolRules List of Boolean rules
#' @param filename Output file path
generateTextReport <- function(ruleStats, boolRules, filename) {
  
  sink(filename)
  
  cat("=============================================================================\n")
  cat("BOOLEAN RULE INFERENCE SUMMARY REPORT\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=============================================================================\n\n")
  
  cat("OVERALL STATISTICS:\n")
  cat("- Total genes with rules:", nrow(ruleStats), "\n")
  cat("- Mean rule quality:", round(mean(ruleStats$bestScore, na.rm = TRUE), 3), "\n")
  cat("- Median rule quality:", round(median(ruleStats$bestScore, na.rm = TRUE), 3), "\n")
  cat("- Rules with quality > 0.8:", sum(ruleStats$bestScore > 0.8, na.rm = TRUE), 
      "(", round(100 * sum(ruleStats$bestScore > 0.8, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n")
  cat("- Rules with quality > 0.95:", sum(ruleStats$bestScore > 0.95, na.rm = TRUE), 
      "(", round(100 * sum(ruleStats$bestScore > 0.95, na.rm = TRUE) / nrow(ruleStats), 1), "%)\n\n")
  
  cat("NETWORK COMPLEXITY:\n")
  cat("- Genes with 1 regulator:", sum(ruleStats$nRegulators == 1, na.rm = TRUE), "\n")
  cat("- Genes with 2 regulators:", sum(ruleStats$nRegulators == 2, na.rm = TRUE), "\n")
  cat("- Genes with 3+ regulators:", sum(ruleStats$nRegulators >= 3, na.rm = TRUE), "\n")
  cat("- Mean regulators per gene:", round(mean(ruleStats$nRegulators, na.rm = TRUE), 2), "\n\n")
  
  if ("methodUsed" %in% colnames(ruleStats)) {
    cat("METHOD USAGE:\n")
    methodCounts <- table(ruleStats$methodUsed, useNA = "no")
    for (method in names(methodCounts)) {
      cat("-", method, ":", methodCounts[method], "genes\n")
    }
    cat("\n")
  }
  
  cat("TOP 10 HIGHEST QUALITY RULES:\n")
  topRules <- ruleStats[order(ruleStats$bestScore, decreasing = TRUE, na.last = TRUE)[1:min(10, nrow(ruleStats))], ]
  for (i in 1:nrow(topRules)) {
    cat(i, ".", topRules$gene[i], "- Score:", round(topRules$bestScore[i], 3), 
        "- Regulators:", topRules$nRegulators[i], "\n")
  }
  
  sink()
  
  message("Text report saved to: ", filename)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x