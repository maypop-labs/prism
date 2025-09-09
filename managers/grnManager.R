# =============================================================================
# grnManager.R (Cleaned Version)
# Functions for creating, visualizing, and saving gene regulatory networks
# Simplified for the new modular GRN pipeline
# =============================================================================

# Helpers
`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Prepare GRN data by filtering expression matrix and setting up SCENIC
#'
#' This function handles the initial data preparation phase for GRN construction,
#' including expression matrix filtering, SCENIC initialization, and gene selection
#' to focus on switch genes and transcription factors.
#'
#' @param cds CellDataSet object from Monocle3 containing expression data
#' @param switchGenes Data frame with switch gene information (from GeneSwitches)
#' @param config Configuration list containing GRN parameters and paths
#' @return List containing filtered expression matrices, SCENIC options, and metadata
#' @export
prepareGrnData <- function(cds, switchGenes, config) {
  if (config$verbose) message("Phase 1: Preparing expression data and SCENIC setup")
  
  # Extract and filter expression matrix
  exprMat <- sparseToDense(assay(cds, "counts"))
  switchGeneNames <- rownames(switchGenes)
  
  # Load TF list from SCENIC
  data(defaultDbNames)
  species <- config$grnScenicSpecies
  baseName <- paste0("motifAnnotations_", species)
  
  # Try to load motif annotations
  candidates <- c(baseName, paste0(baseName, "_v10"), paste0(baseName, "_v9"), "motifAnnotations")
  loaded <- FALSE
  for (ds in candidates) {
    if (ds %in% data(package = "RcisTarget")$results[, "Item"]) {
      data(list = ds, package = "RcisTarget", envir = environment())
      assign(baseName, get(ds), envir = globalenv())
      loaded <- TRUE
      break
    }
  }
  if (!loaded) stop("No motif annotation found for '", species, "'")
  
  # Initialize SCENIC options
  scenicOptions <- initializeScenic(
    org = config$grnScenicSpecies,
    dbDir = config$rcisTargetPath, 
    dbs = config$grnScenicDBs,
    nCores = config$cores
  )
  
  # Filter genes
  allTFs <- getDbTfs(scenicOptions)
  keepGenes <- rowSums(exprMat > 1) >= 10
  exprMat <- exprMat[keepGenes, ]
  
  # Focus on switch genes + TFs
  unionGenes <- union(switchGeneNames, allTFs)
  exprMat <- exprMat[intersect(rownames(exprMat), unionGenes), ]
  
  # SCENIC gene filtering
  filteredGenes <- geneFiltering(exprMat, scenicOptions)
  exprMatLog <- log2(exprMat[filteredGenes, ] + 1)
  
  if (config$verbose) {
    message("Expression matrix: ", nrow(exprMatLog), " genes x ", ncol(exprMatLog), " cells")
    message("Switch genes in matrix: ", sum(switchGeneNames %in% rownames(exprMatLog)))
  }
  
  return(list(
    exprMat = exprMat,
    exprMatLog = exprMatLog,
    scenicOptions = scenicOptions,
    filteredGenes = filteredGenes,
    switchGenes = switchGenes
  ))
}

#' Add GeneSwitches temporal information to regulatory edges
#'
#' Integrates temporal information from GeneSwitches analysis with SCENIC edges,
#' including switch timing, direction, and quality scores. Also computes temporal
#' consistency metrics to identify edges where TFs switch before their targets.
#'
#' @param edges Data frame with TF and Target columns (from SCENIC)
#' @param switchGenes Data frame with switch gene analysis results
#' @param verbose Logical, whether to print progress information
#' @return Enhanced edges data frame with temporal information columns
#' @export
addGeneSwitchesInfo <- function(edges, switchGenes, verbose = FALSE) {
  if (verbose) message("  Adding GeneSwitches temporal information...")
  
  # Add switch timing for targets
  edges$targetSwitchTime <- NA_real_
  edges$targetDirection <- NA_character_
  edges$targetPseudoR2 <- NA_real_
  
  matchedTargets <- match(edges$Target, rownames(switchGenes))
  validMatches <- !is.na(matchedTargets)
  
  if (any(validMatches)) {
    edges$targetSwitchTime[validMatches] <- switchGenes$switch_at_time[matchedTargets[validMatches]]
    edges$targetDirection[validMatches] <- switchGenes$direction[matchedTargets[validMatches]] 
    edges$targetPseudoR2[validMatches] <- switchGenes$pseudoR2s[matchedTargets[validMatches]]
  }
  
  # Add switch timing for TFs (if available)
  edges$tfSwitchTime <- NA_real_
  edges$tfDirection <- NA_character_
  
  matchedTFs <- match(edges$TF, rownames(switchGenes))
  validTFMatches <- !is.na(matchedTFs)
  
  if (any(validTFMatches)) {
    edges$tfSwitchTime[validTFMatches] <- switchGenes$switch_at_time[matchedTFs[validTFMatches]]
    edges$tfDirection[validTFMatches] <- switchGenes$direction[matchedTFs[validTFMatches]]
  }
  
  # Calculate temporal consistency (TF should switch before target)
  edges$temporallyConsistent <- NA
  bothHaveTiming <- !is.na(edges$tfSwitchTime) & !is.na(edges$targetSwitchTime)
  if (any(bothHaveTiming)) {
    edges$temporallyConsistent[bothHaveTiming] <- 
      edges$tfSwitchTime[bothHaveTiming] <= edges$targetSwitchTime[bothHaveTiming]
  }
  
  if (verbose) {
    nTargetsWithTiming <- sum(!is.na(edges$targetSwitchTime))
    nTFsWithTiming <- sum(!is.na(edges$tfSwitchTime))
    nTemporallyConsistent <- sum(edges$temporallyConsistent, na.rm = TRUE)
    message("    Targets with switch timing: ", nTargetsWithTiming, "/", nrow(edges))
    message("    TFs with switch timing: ", nTFsWithTiming, "/", nrow(edges))
    message("    Temporally consistent edges: ", nTemporallyConsistent)
  }
  
  return(edges)
}

#' Extract SCENIC metadata from regulonTargetsInfo
#' 
#' Simplified version that assumes standard SCENIC output structure
#' and focuses on essential motif information.
#' 
#' @param scenicOptions SCENIC options object (used for file paths)
#' @param edges Data frame with TF and Target columns
#' @param verbose Logical, whether to print progress messages
#' @return Enhanced edges data frame with motif metadata
#' @export
extractScenicMetadata <- function(scenicOptions, edges, verbose = TRUE) {
  if (verbose) message("Extracting SCENIC motif metadata...")

  # Direct file loading
  rtiFile <- "int/2.5_regulonTargetsInfo.Rds"
  if (!file.exists(rtiFile)) {
    warning("regulonTargetsInfo file not found, adding empty motif columns")
    edges$hasMotif <- FALSE
    edges$motifConfidence <- 0
    edges$NES <- NA_real_
    edges$nMotifs <- NA_integer_
    edges$bestMotif <- NA_character_
    return(edges)
  }
  
  rti <- readRDS(rtiFile)
  rti <- as.data.frame(rti)  # Convert data.table to data.frame to avoid merge conflicts

  # Rename gene column to Target for consistency BEFORE selecting columns
  if ("gene" %in% names(rti)) {
    names(rti)[names(rti) == "gene"] <- "Target"
  }

  # Select essential columns
  essentialCols <- c("TF", "Target", "NES", "nMotifs", "highConfAnnot", "bestMotif")
  availableCols <- intersect(essentialCols, names(rti))
  
  if (!"TF" %in% availableCols || !"Target" %in% availableCols) {
    warning("Required TF/Target columns missing from regulonTargetsInfo")
    edges$hasMotif <- FALSE
    edges$motifConfidence <- 0
    return(edges)
  }
  
  # Keep highest NES per TF-Target pair if duplicates exist
  if ("NES" %in% names(rti)) {
    rti <- rti[order(-ifelse(is.finite(rti$NES), rti$NES, -Inf)), ]
  }
  
  # Remove duplicates
  rti <- rti[!duplicated(rti[, c("TF", "Target")]), ]
  
  # Merge with edges
  result <- merge(edges, rti[availableCols], by = c("TF", "Target"), all.x = TRUE)
  
  # After merge, ensure numeric columns are actually numeric
  if ("NES" %in% names(result)) {
    result$NES <- as.numeric(result$NES)
  }
  if ("nMotifs" %in% names(result)) {
    result$nMotifs <- as.numeric(result$nMotifs)
  }
  
  # Create derived motif columns
  result$hasMotif <- FALSE
  result$motifConfidence <- 0
  
  if ("NES" %in% names(result)) {
    result$hasMotif <- is.finite(result$NES) & result$NES > 2.0
    result$motifConfidence <- pmin(1, pmax(0, (result$NES - 2.0) / 3.0))
  }
  
  # Fallback for missing NES
  if ("nMotifs" %in% names(result)) {
    noNES <- is.na(result$NES) | !is.finite(result$NES)
    result$hasMotif[noNES] <- !is.na(result$nMotifs[noNES]) & result$nMotifs[noNES] > 0
    result$motifConfidence[noNES] <- ifelse(result$hasMotif[noNES], 0.6, 0)
  }
  
  if (verbose) {
    nWithMotif <- sum(result$hasMotif, na.rm = TRUE)
    message("  Edges with motif support: ", nWithMotif, "/", nrow(result))
  }
  
  return(result)
}

#' Assign edge signs using motif-aware approach
#'
#' Prioritizes correlation over motif when they conflict, as correlation
#' reflects observed data while motif conflicts may indicate hidden regulators.
#'
#' @param edges Data frame with correlation and motif information  
#' @param config Configuration object with thresholds
#' @param verbose Logical, whether to print summary
#' @return Edges with regType column added
#' @export
assignMotifAwareSigns <- function(edges, config, verbose = FALSE) {
  
  if (verbose) message("Assigning regulatory signs...")
  
  # Remove edges with missing correlation data
  validCorr <- !is.na(edges$corr) & is.finite(edges$corr)
  if (!all(validCorr)) {
    nRemoved <- sum(!validCorr)
    edges <- edges[validCorr, ]
    if (verbose) message("  Removed ", nRemoved, " edges with invalid correlations")
  }
  
  # Ensure motifConfidence exists and handle NAs
  if (!"motifConfidence" %in% colnames(edges)) {
    edges$motifConfidence <- 0
  }
  edges$motifConfidence[is.na(edges$motifConfidence)] <- 0
  
  edges$regType <- case_when(
    # Strong correlations override motifs (data-driven approach)
    edges$corr > config$grnPositiveThreshold ~ "Activation",
    edges$corr < -config$grnNegativeThreshold ~ "Inhibition",
    
    # High motif confidence for weak correlations
    edges$motifConfidence > config$grnHighMotifThreshold & 
      edges$corr > config$grnMotifCorrAlignThreshold ~ "Activation",
    
    edges$motifConfidence > config$grnHighMotifThreshold & 
      edges$corr < -config$grnMotifCorrAlignThreshold ~ "Inhibition",
    
    # Medium motif confidence as tiebreaker
    edges$motifConfidence > config$grnMediumMotifThreshold & 
      edges$corr >= 0 ~ "Activation",
    
    edges$motifConfidence > config$grnMediumMotifThreshold & 
      edges$corr < 0 ~ "Inhibition",
    
    # Default: follow correlation direction
    TRUE ~ ifelse(edges$corr >= 0, "Activation", "Inhibition")
  )
  
  if (verbose) {
    signCounts <- table(edges$regType)
    message("  Activation edges: ", signCounts["Activation"] %||% 0)
    message("  Inhibition edges: ", signCounts["Inhibition"] %||% 0)
  }
  
  return(edges)
}

#' Compute composite ranking for edge importance
#'
#' Simplified version with automatic weight normalization.
#'
#' @param edges Data frame with correlation, GENIE3, and motif scores
#' @param config Configuration object with weights
#' @param verbose Logical, whether to print progress  
#' @return Edges with composite scoring columns
#' @export
computeCompositeRanking <- function(edges, config, verbose = FALSE) {
  if (verbose) message("Computing composite edge ranking...")
  
  # Normalize weights
  weights <- c(config$grnCorrWeight, config$grnGenie3Weight, config$grnMotifWeight)
  weights <- weights / sum(weights)
  
  # Normalize individual scores to [0,1] with division by zero guards
  maxCorr <- max(abs(edges$corr), na.rm = TRUE)
  maxGenie3 <- max(edges$genie3Importance, na.rm = TRUE)
  
  # Guard against division by zero
  if (maxCorr <= 0 || !is.finite(maxCorr)) {
    edges$corrScore <- 0.1  # Small default when no correlations exist
    if (verbose) message("  Warning: No valid correlations found, using default correlation scores")
  } else {
    edges$corrScore <- abs(edges$corr) / maxCorr
  }
  
  if (maxGenie3 <= 0 || !is.finite(maxGenie3)) {
    edges$genie3Score <- 0.1  # Small default when no GENIE3 scores exist
    if (verbose) message("  Warning: No valid GENIE3 scores found, using default GENIE3 scores")
  } else {
    edges$genie3Score <- edges$genie3Importance / maxGenie3
  }
  
  # Handle motif scores (already bounded)
  edges$motifScore <- pmax(0, pmin(1, edges$motifConfidence))
  edges$motifScore[is.na(edges$motifScore)] <- 0
  
  # Composite score
  edges$compositeScore <- weights[1] * edges$corrScore + 
                         weights[2] * edges$genie3Score + 
                         weights[3] * edges$motifScore
  
  # Rank by target gene
  edges <- edges %>%
    group_by(Target) %>%
    mutate(regulatorRank = dense_rank(dplyr::desc(compositeScore))) %>%
    ungroup()
  
  if (verbose) {
    avgScore <- mean(edges$compositeScore, na.rm = TRUE)
    message("  Average composite score: ", round(avgScore, 3))
  }
  
  return(edges)
}

#' Convert sparse matrix to dense for SCENIC compatibility
#' 
#' @param exprMat Expression matrix (potentially sparse)
#' @return Dense matrix
#' @export
sparseToDense <- function(exprMat) {
  if (methods::is(exprMat, "dgCMatrix") || methods::is(exprMat, "Matrix")) {
    exprMat <- as.matrix(exprMat)
  }
  exprMat
}

#' Create and save a GRN plot with enhanced metadata
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param plotPath character, full path where plot should be saved
#' @param title character, plot title (default: "Gene Regulatory Network")
#' @param config list, configuration object containing plot parameters
#' @param edges data.frame, edge metadata for enhanced visualization
#' @param verbose logical, whether to print status messages
#' @return ggplot object (invisibly)
#' @export
createAndSaveGrnPlot <- function(graph, plotPath, title = "Gene Regulatory Network", 
                                 config, edges = NULL, verbose = FALSE) {
  
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (igraph::vcount(graph) == 0) {
    warning("Graph has no vertices - creating empty plot")
    return(invisible(NULL))
  }
  
  if (verbose) {
    message("Creating enhanced GRN plot with ", igraph::vcount(graph), " nodes and ", igraph::ecount(graph), " edges")
  }
  
  # Add edge attributes if available
  if (!is.null(edges)) {
    edgeList <- igraph::as_data_frame(graph, what = "edges")
    enhancedEdges <- edgeList %>%
      left_join(edges %>% select(from = TF, to = Target, regType, motifConfidence), 
                by = c("from", "to"))
    
    # Add to graph
    if ("regType" %in% colnames(enhancedEdges)) {
      igraph::E(graph)$regType <- enhancedEdges$regType
    }
    if ("motifConfidence" %in% colnames(enhancedEdges)) {
      igraph::E(graph)$motifConfidence <- enhancedEdges$motifConfidence
    }
  }
  
  # Create enhanced plot with curved edges and improved styling
  grnPlot <- ggraph::ggraph(graph, layout = config$networkLayout) +
    {
      if (igraph::ecount(graph) > 0) {
        if ("regType" %in% igraph::edge_attr_names(graph)) {
          ggraph::geom_edge_arc(aes(colour = regType), 
                               start_cap = ggraph::circle(0, 'mm'),
                               end_cap = ggraph::circle(6, 'mm'),
                               alpha = 1.0,
                               width = 1.2,
                               curvature = 0.2,
                               arrow = grid::arrow(length = grid::unit(4, "mm"), 
                                                  type = "closed"))
        } else {
          ggraph::geom_edge_arc(alpha = 1.0,
                               start_cap = ggraph::circle(0, 'mm'),
                               end_cap = ggraph::circle(6, 'mm'),
                               width = 1.2,
                               curvature = 0.2,
                               arrow = grid::arrow(length = grid::unit(4, "mm")))
        }
      } else {
        ggraph::geom_edge_arc(alpha = 1.0)
      }
    } +
    ggraph::geom_node_point(size = 8, 
                           color = "black",
                           stroke = 1.5,
                           shape = 21,
                           fill = "lightblue") +
    ggraph::geom_node_text(aes(label = name), 
                          size = 3.5, 
                          repel = TRUE,
                          point.padding = 1.2,
                          box.padding = 0.8,
                          force = 5,
                          max.overlaps = Inf) +
    ggraph::scale_edge_colour_manual(values = c("Activation" = "darkgreen", 
                                               "Inhibition" = "darkred"),
                                     guide = "none") +
    ggplot2::ggtitle(title) +
    ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    ggraph::theme_graph()
  
  # Save plot
  if (verbose) message("Saving enhanced GRN plot to: ", plotPath)
  
  ggplot2::ggsave(
    filename = plotPath,
    plot = grnPlot,
    width = config$figWidth,
    height = config$figHeight,
    dpi = config$figDPI,
    units = "in"
  )
  
  invisible(grnPlot)
}

#' Create and save a GraphML file with enhanced metadata
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param graphmlPath character, full path where GraphML file should be saved
#' @param edges data.frame, edge metadata to include in GraphML
#' @param verbose logical, whether to print status messages
#' @return NULL (invisibly)
#' @export
createAndSaveGraphml <- function(graph, graphmlPath, edges = NULL, verbose = FALSE) {
  
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (igraph::vcount(graph) == 0 || igraph::ecount(graph) == 0) {
    stop("Cannot save GraphML: graph has no vertices or edges")
  }
  
  if (verbose) message("Saving GraphML file to: ", graphmlPath)
  
  # Add edge attributes if available
  if (!is.null(edges)) {
    edgeList <- igraph::as_data_frame(graph, what = "edges")
    enhancedEdges <- edgeList %>%
      left_join(edges %>% 
                select(from = TF, to = Target, 
                       any_of(c("genie3Importance", "regType", "motifConfidence", "corr"))), 
                by = c("from", "to"))
    
    # Add attributes to graph
    for (col in c("genie3Importance", "regType", "motifConfidence", "corr")) {
      if (col %in% colnames(enhancedEdges)) {
        if (is.numeric(enhancedEdges[[col]])) {
          enhancedEdges[[col]][is.na(enhancedEdges[[col]])] <- 0
        } else {
          enhancedEdges[[col]][is.na(enhancedEdges[[col]])] <- "Unknown"
        }
        igraph::E(graph)[[col]] <- enhancedEdges[[col]]
      }
    }
  }
  
  # Ensure directory exists
  dir.create(dirname(graphmlPath), recursive = TRUE, showWarnings = FALSE)
  
  # Save GraphML
  tryCatch({
    igraph::write_graph(graph, graphmlPath, format = "graphml")
    if (verbose) {
      message("Successfully saved GraphML with ", igraph::vcount(graph), 
              " nodes and ", igraph::ecount(graph), " edges")
    }
  }, error = function(e) {
    stop("Failed to save GraphML file: ", e$message)
  })
  
  invisible(NULL)
}

#' Save complete enhanced GRN output set
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param edges data.frame with enhanced edge metadata
#' @param rdsPath character, path for RDS file
#' @param plotPath character, path for plot file
#' @param graphmlPath character, path for GraphML file
#' @param edgesPath character, path for edges CSV file
#' @param title character, plot title
#' @param config list, configuration object
#' @param verbose logical, whether to print status messages
#' @return list with paths to saved files (invisibly)
#' @export
saveEnhancedGrnOutputSet <- function(graph, edges, rdsPath, plotPath, graphmlPath, 
                                     edgesPath, title, config, verbose = FALSE) {
  
  if (verbose) message("Saving complete GRN output set...")
  
  # Ensure directories exist
  for (path in c(rdsPath, plotPath, graphmlPath, edgesPath)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save RDS
  if (verbose) message("  Saving GRN object...")
  saveRDS(graph, file = rdsPath)
  
  # Save edges CSV with all computed columns
  if (verbose) message("  Saving comprehensive edges table...")
  
  # Select all relevant columns for export
  exportCols <- c("TF", "Target", "regType", "corr", "genie3Importance", 
                  "motifConfidence", "corrScore", "genie3Score", "motifScore", 
                  "compositeScore", "regulatorRank")
  
  # Add temporal columns if available
  temporalCols <- c("targetSwitchTime", "targetDirection", "targetPseudoR2", 
                    "tfSwitchTime", "tfDirection", "temporallyConsistent")
  exportCols <- c(exportCols, intersect(temporalCols, names(edges)))
  
  # Add SCENIC metadata if available
  scenicCols <- c("NES", "nMotifs", "bestMotif", "hasMotif")
  exportCols <- c(exportCols, intersect(scenicCols, names(edges)))
  
  # Export only available columns
  availableCols <- intersect(exportCols, names(edges))
  if (length(availableCols) < length(exportCols)) {
    missingCols <- setdiff(exportCols, availableCols)
    if (verbose) message("  Note: Missing columns in export: ", paste(missingCols, collapse = ", "))
  }
  
  write.table(edges[availableCols], file = edgesPath, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save plot
  if (verbose) message("  Creating plot...")
  createAndSaveGrnPlot(graph, plotPath, title, config, edges, verbose = FALSE)
  
  # Save GraphML
  if (verbose) message("  Saving GraphML...")
  createAndSaveGraphml(graph, graphmlPath, edges, verbose = FALSE)
  
  savedFiles <- list(
    rds = rdsPath,
    edges = edgesPath,
    plot = plotPath,
    graphml = graphmlPath
  )
  
  if (verbose) {
    message("Completed saving GRN output set with ", length(availableCols), " edge attributes")
  }
  
  invisible(savedFiles)
}

# ----------------------------------------------------------------------
# buildEnhancedGrn
# Constructs enhanced GRN by integrating SCENIC with GeneSwitches
# ----------------------------------------------------------------------

#' Build enhanced gene regulatory network
#'
#' Integrates SCENIC regulons with GeneSwitches temporal information to create
#' an enhanced gene regulatory network. Performs correlation analysis, motif-aware
#' sign assignment, composite ranking, and regulator capping for downstream
#' Boolean network construction.
#'
#' @param scenicOptions SCENIC options object containing pipeline results
#' @param switchGenes Data frame with switch gene analysis results
#' @param exprMat Expression matrix for correlation calculations
#' @param config Configuration object with GRN parameters
#' @return Data frame with enhanced edge metadata
#' @export
buildEnhancedGrn <- function(scenicOptions, switchGenes, exprMat, config) {
  if (config$verbose) message("Phase 3: Building enhanced GRN with GeneSwitches integration")
  
  # Load SCENIC outputs with error handling
  tryCatch({
    regulons <- loadInt(scenicOptions, "regulons")
  }, error = function(e) {
    stop("Failed to load SCENIC regulons: ", e$message, 
         "\nCheck that SCENIC pipeline completed successfully.")
  })
  
  if (length(regulons) == 0) {
    stop("No regulons found - SCENIC pipeline may have failed or no significant regulons detected")
  }
  
  # Extract edges from regulons, filter to switch genes  
  scenicEdges <- purrr::map_dfr(names(regulons), function(tf) {
    targets <- regulons[[tf]]
    if (length(targets) > 0) {
      data.frame(TF = tf, Target = targets, stringsAsFactors = FALSE)
    }
  })
  
  if (nrow(scenicEdges) == 0) {
    stop("No edges extracted from regulons - check SCENIC output integrity")
  }
  
  switchGeneNames <- rownames(switchGenes)
  scenicEdges <- scenicEdges %>% filter(Target %in% switchGeneNames)
  
  if (nrow(scenicEdges) == 0) {
    stop("No edges connect to switch genes - check gene name matching between SCENIC and GeneSwitches outputs")
  }
  
  if (config$verbose) {
    message("Initial SCENIC edges to switch genes: ", nrow(scenicEdges))
  }
  
  # Add SCENIC metadata (motif confidence, NES, etc.)
  scenicEdges <- extractScenicMetadata(scenicOptions, scenicEdges, config$verbose)
  
  # Add correlations with improved error handling
  validGenes <- intersect(rownames(exprMat), unique(c(scenicEdges$TF, scenicEdges$Target)))
  if (length(validGenes) < 2) {
    stop("Insufficient valid genes for correlation analysis - check expression matrix and gene names")
  }
  
  exprMat <- exprMat[validGenes, ]
  scenicEdges <- scenicEdges %>% filter(TF %in% validGenes & Target %in% validGenes)
  
  getCorrelation <- function(tf, tg) {
    tryCatch({
      cor(exprMat[tf, ], exprMat[tg, ], method = "spearman", use = "complete.obs")
    }, error = function(e) {
      warning("Correlation failed for ", tf, " -> ", tg, ": ", e$message)
      return(0)
    })
  }
  scenicEdges$corr <- mapply(getCorrelation, scenicEdges$TF, scenicEdges$Target)
  
  # Remove edges with failed correlations
  validCorr <- is.finite(scenicEdges$corr)
  if (!all(validCorr)) {
    nRemoved <- sum(!validCorr)
    scenicEdges <- scenicEdges[validCorr, ]
    if (config$verbose) message("Removed ", nRemoved, " edges with invalid correlations")
  }
  
  # Assign regulatory signs using motif-aware approach
  scenicEdges <- assignMotifAwareSigns(scenicEdges, config, config$verbose)
  
  # Add GeneSwitches temporal information
  scenicEdges <- addGeneSwitchesInfo(scenicEdges, switchGenes, config$verbose)
  
  # Compute composite ranking
  scenicEdges <- computeCompositeRanking(scenicEdges, config, config$verbose)
  
  # Cap regulators per target gene to improve Boolean rule learning
  maxRegulatorsPerTarget <- config$boolMaxRegulators # From config for Boolean learning
  scenicEdges <- scenicEdges %>%
    group_by(Target) %>%
    arrange(desc(compositeScore)) %>%
    slice_head(n = maxRegulatorsPerTarget) %>%
    ungroup()
  
  if (config$verbose) {
    message("Applied regulator capping: ", nrow(scenicEdges), " edges retained (max ", 
            maxRegulatorsPerTarget, " regulators per target)")
  }
  
  return(scenicEdges)
}

# ----------------------------------------------------------------------  
# filterAndFinalizeGrn
# Applies filtering criteria and creates final igraph object
# ----------------------------------------------------------------------

#' Filter and finalize gene regulatory network
#'
#' Applies comprehensive filtering using correlation, motif, temporal, and
#' quality criteria. Creates final igraph object and optionally removes
#' isolated nodes. Provides detailed filtering statistics.
#'
#' @param edges Data frame with enhanced edge metadata
#' @param config Configuration object with filtering thresholds
#' @return List containing final igraph object and filtered edges
#' @export
filterAndFinalizeGrn <- function(edges, config) {
  if (config$verbose) message("Phase 4: Filtering and finalizing GRN")
  
  # Remove edges with critical NA values first
  criticalCols <- c("corr", "TF", "Target")
  validRows <- complete.cases(edges[criticalCols])
  if (!all(validRows)) {
    nRemoved <- sum(!validRows)
    edges <- edges[validRows, ]
    if (config$verbose) message("Removed ", nRemoved, " edges with missing critical data")
  }
  
  if (nrow(edges) == 0) {
    stop("No edges remain after NA filtering - check data quality")
  }
  
  # Apply enhanced filtering with temporal and motif priors
  keepEdges <- with(edges, {
    # Basic correlation filters (handle NAs safely)
    strongCorr <- (!is.na(corr)) & ((corr > config$grnPositiveThreshold) | (corr < -config$grnNegativeThreshold))
    
    # Motif-based priors (use motifConfidence if available)
    motifPrior <- FALSE
    if ("motifConfidence" %in% names(edges)) {
      motifPrior <- (!is.na(motifConfidence)) & (motifConfidence > config$grnPriorThreshold)
    }
    
    # Temporal consistency bonus (handle NAs)
    temporalPrior <- FALSE  
    if ("temporallyConsistent" %in% names(edges)) {
      temporalPrior <- (!is.na(temporallyConsistent)) & (temporallyConsistent == TRUE)
    }
    
    # High-quality switch genes get priority (handle NAs)
    highQualityTarget <- FALSE
    if ("targetPseudoR2" %in% names(edges)) {
      highQualityTarget <- (!is.na(targetPseudoR2)) & (targetPseudoR2 > config$grnSwitchQualityThreshold)
    }
    
    # Keep edges that meet any of these criteria
    strongCorr | motifPrior | temporalPrior | highQualityTarget
  })
  
  if (!any(keepEdges, na.rm = TRUE)) {
    stop("No edges passed filtering - consider relaxing thresholds in config")
  }
  
  edges <- edges[keepEdges, ]
  
  if (config$verbose) {
    message("Retained ", nrow(edges), " edges after enhanced filtering")
    
    # Print filtering statistics
    nStrongCorr <- sum((!is.na(edges$corr)) & ((edges$corr > config$grnPositiveThreshold) | (edges$corr < -config$grnNegativeThreshold)))
    message("  - Strong correlation: ", nStrongCorr, " edges")
    
    if ("motifConfidence" %in% names(edges)) {
      nMotif <- sum((!is.na(edges$motifConfidence)) & (edges$motifConfidence > config$grnPriorThreshold))
      message("  - Motif support: ", nMotif, " edges")
    }
  }
  
  # Build igraph with error handling
  tryCatch({
    g <- graph_from_data_frame(edges, directed = TRUE)
  }, error = function(e) {
    stop("Failed to create igraph object: ", e$message, "\nCheck edge data format")
  })
  
  # Optional: remove isolated nodes
  if (config$grnRemoveIsolatedNodes) {
    isolatedNodes <- names(which(degree(g, mode = "all") == 0))
    if (length(isolatedNodes) > 0) {
      g <- delete_vertices(g, isolatedNodes) 
      edges <- edges %>% filter(TF %in% V(g)$name & Target %in% V(g)$name)
      if (config$verbose) message("Removed ", length(isolatedNodes), " isolated nodes")
    }
  }
  
  if (config$verbose) {
    message("Final GRN: ", vcount(g), " nodes, ", ecount(g), " edges")
  }
  
  return(list(graph = g, edges = edges))
}

# ----------------------------------------------------------------------
# runGenie3Pipeline  
# Orchestrates GENIE3 execution with checkpoint recovery
# ----------------------------------------------------------------------

#' Run GENIE3 pipeline with checkpoint support
#'
#' Manages the GENIE3 inference pipeline including correlation matrix computation
#' and GENIE3 weight matrix calculation. Supports checkpoint recovery to avoid
#' re-running expensive computations.
#'
#' @param prepData List containing prepared expression data and SCENIC options
#' @param config Configuration object with pipeline settings
#' @return Updated SCENIC options object
#' @export
runGenie3Pipeline <- function(prepData, config) {
  if (config$verbose) message("Phase 2: Running GENIE3 pipeline")
  
  # Check if GENIE3 already completed
  genie3File <- "int/1.4_GENIE3_linkList.Rds"
  if (file.exists(genie3File) && !config$grnFromScratch) {
    if (config$verbose) message("Found existing GENIE3 results, loading...")
    return(prepData$scenicOptions)
  }
  
  if (config$verbose) message("Running GENIE3 from scratch...")
  
  # Step 2.1: Correlation matrix (checkpoint)
  corrFile <- "int/1.2_corrMat.Rds" 
  if (!file.exists(corrFile) || config$grnFromScratch) {
    if (config$verbose) message("  Computing correlation matrix...")
    runCorrelation(prepData$exprMatLog, prepData$scenicOptions)
  }
  
  # Step 2.2: GENIE3 inference (checkpoint) 
  genieFiles <- list.files("int", pattern = "GENIE3_weightMatrix_part_", full.names = TRUE)
  if (length(genieFiles) == 0 || config$grnFromScratch) {
    if (config$verbose) message("  Running GENIE3 (this will take several hours)...")
    runGenie3(prepData$exprMatLog, prepData$scenicOptions)
  } else {
    if (config$verbose) message("  Found existing GENIE3 weight matrices")
  }
  
  return(prepData$scenicOptions)
}

# ----------------------------------------------------------------------
# runScenicModules
# Handles SCENIC module and regulon creation phases
# ----------------------------------------------------------------------

#' Run SCENIC modules and regulon creation
#'
#' Manages the SCENIC pipeline phases for creating co-expression modules
#' and regulons with motif enrichment. Supports checkpoint recovery.
#'
#' @param scenicOptions SCENIC options object from GENIE3 phase
#' @param config Configuration object with pipeline settings
#' @return Updated SCENIC options object
#' @export
runScenicModules <- function(scenicOptions, config) {
  if (config$verbose) message("Phase 2b: Creating TF modules and regulons")
  
  # Check if modules already exist
  modulesFile <- "int/1.6_tfModules_asDF.Rds"
  if (file.exists(modulesFile) && !config$grnFromScratch) {
    if (config$verbose) message("Found existing TF modules")
  } else {
    if (config$verbose) message("  Creating co-expression modules...")
    scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  }
  
  # Check if regulons already exist  
  regulonsFile <- "int/2.6_regulons_asGeneSet.Rds"
  if (file.exists(regulonsFile) && !config$grnFromScratch) {
    if (config$verbose) message("Found existing regulons")
  } else {
    if (config$verbose) message("  Creating regulons with motif enrichment...")
    scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, onlyPositiveCorr = config$grnOnlyPositiveCorr)
  }
  
  return(scenicOptions)
}

# ----------------------------------------------------------------------
# runScenicScoring
# Manages SCENIC AUCell scoring phase
# ----------------------------------------------------------------------

#' Run SCENIC regulon scoring
#'
#' Manages the SCENIC AUCell scoring phase to compute regulon activity
#' scores across cells. Required for complete SCENIC metadata extraction.
#'
#' @param scenicOptions SCENIC options object from regulon creation phase
#' @param exprMatLog Log-transformed expression matrix
#' @param config Configuration object with pipeline settings
#' @return Updated SCENIC options object
#' @export
runScenicScoring <- function(scenicOptions, exprMatLog, config) {
  if (config$verbose) message("Phase 2c: Scoring regulons on cells")
  
  # Check if AUC already exists
  aucFile <- "int/3.4_regulonAUC.Rds" 
  if (file.exists(aucFile) && !config$grnFromScratch) {
    if (config$verbose) message("Found existing regulon AUC scores")
  } else {
    if (config$verbose) message("  Computing AUCell scores...")
    scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMatLog)
  }
  
  return(scenicOptions)
}

# ----------------------------------------------------------------------
# saveGrnOutputs
# Coordinates saving of all GRN outputs and cleanup
# ----------------------------------------------------------------------

#' Save GRN outputs and perform cleanup
#'
#' Coordinates saving of all GRN outputs including RDS, plots, GraphML,
#' and TSV files. Optionally preserves SCENIC intermediate files for
#' debugging or future analysis.
#'
#' @param grn List containing final igraph object and edges data frame
#' @param ptPaths Trajectory-specific paths from getTrajectoryFilePaths
#' @param config Configuration object with save settings
#' @return GRN object (invisibly)
#' @export
saveGrnOutputs <- function(grn, ptPaths, config) {
  if (config$verbose) message("Phase 5: Saving outputs and cleaning up")
  
  if (!config$saveResults) {
    if (config$verbose) message("saveResults = FALSE, skipping file output")
    return(invisible(grn))
  }
  
  # Use trajectory-specific output paths
  rdsPath     <- ptPaths$grnEdges
  plotPath    <- ptPaths$grnPlot
  graphmlPath <- ptPaths$grnGraphml
  edgesPath   <- ptPaths$grnEdgesTsv
  
  # Save complete output set
  saveEnhancedGrnOutputSet(
    graph       = grn$graph,
    edges       = grn$edges, 
    rdsPath     = rdsPath,
    plotPath    = plotPath,
    graphmlPath = graphmlPath,
    edgesPath   = edgesPath,
    title       = "Gene Regulatory Network (GRN)",
    config      = config,
    verbose     = config$verbose
  )
  
  # Optional: preserve key SCENIC intermediate files
  if (config$grnPreserveIntermediates) {
    tryCatch({
      if (dir.exists("int")) {
        file.copy("int", ptPaths$scenic, recursive = TRUE)
      }
      if (dir.exists("output")) {
        file.copy("output", ptPaths$scenic, recursive = TRUE)
      }
      if (config$verbose) message("Preserved SCENIC files to: ", ptPaths$scenic)
      
      # Clean up original folders after successful preservation
      if (config$grnCleanupAfterPreservation) {
        if (dir.exists("int")) {
          unlink("int", recursive = TRUE)
          if (config$verbose) message("Cleaned up original int/ folder")
        }
        if (dir.exists("output")) {
          unlink("output", recursive = TRUE)
          if (config$verbose) message("Cleaned up original output/ folder")
        }
      }
    }, error = function(e) {
      warning("Could not preserve SCENIC files: ", e$message)
    })
  }
  
  invisible(grn)
}
