# =============================================================================
# grnManager.R
# Functions for creating, visualizing, and saving gene regulatory networks
# with improved integration of SCENIC/GENIE3 biological priors
# =============================================================================

#' Preflight validation of SCENIC data before long runs
#' 
#' @param scenicOptions SCENIC options object
#' @param scenicEdges data.frame with TF, Target columns
#' @return TRUE if validation passes, stops with error if not
preflightPrism <- function(scenicOptions, scenicEdges) {
  
  cat("=== PREFLIGHT VALIDATION ===\n")
  
  # Check TF database
  allTFs <- getDbTfs(scenicOptions)
  cat("TF database size:", length(allTFs), "\n")
  if (length(allTFs) < 100) {
    stop("PREFLIGHT FAILED: Too few TFs in database (", length(allTFs), ")")
  }
  
  # Check regulonTargetsInfo
  tryCatch({
    ri <- loadInt(scenicOptions, "regulonTargetsInfo")
    cat("RegulonTargetsInfo rows:", nrow(ri), "\n")
    cat("Unique TFs in RI:", length(unique(ri$TF)), "\n")
    
    # Check TF recognition
    recog <- mean(unique(ri$TF) %in% allTFs)
    cat(sprintf("TF recognition vs DB: %.1f%%\n", 100*recog))
    if (recog < 0.8) {
      stop("PREFLIGHT FAILED: Low TF recognition rate (", round(recog*100,1), "%). Check symbol sets/organism/annotation.")
    }
    
    # Check required columns (handle gene/target variation)
    if ("gene" %in% colnames(ri) && !"target" %in% colnames(ri)) {
      cat("Note: regulonTargetsInfo uses 'gene' column (will be renamed to 'target')\n")
      requiredCols <- c("TF", "gene")
    } else {
      requiredCols <- c("TF", "target") 
    }
    
    if (!all(requiredCols %in% colnames(ri))) {
      stop("PREFLIGHT FAILED: regulonTargetsInfo missing required columns: ", paste(setdiff(requiredCols, colnames(ri)), collapse=", "))
    }
    
  }, error = function(e) {
    stop("PREFLIGHT FAILED: Cannot load regulonTargetsInfo: ", e$message)
  })
  
  # Dry-run the extractor on small sample
  if (nrow(scenicEdges) > 1000) {
    edgesSmall <- scenicEdges %>% slice_head(n = 1000)
  } else {
    edgesSmall <- scenicEdges
  }
  
  cat("Testing extractor on", nrow(edgesSmall), "edges...\n")
  suppressMessages({
    tryCatch({
      e <- extractScenicMetadata(scenicOptions, scenicEdges = edgesSmall, verbose = FALSE)
      cat("✓ Preflight edges processed:", nrow(e), "\n")
      
      # Belt-and-suspenders: ensure critical columns exist and are valid
      stopifnot("motifNES" %in% colnames(e))
      stopifnot("motifConfidence" %in% colnames(e))
      stopifnot(is.numeric(e$motifNES), is.numeric(e$motifConfidence))
      stopifnot(nrow(e) > 0)
      
    }, error = function(err) {
      stop("PREFLIGHT FAILED: Extractor test failed: ", err$message)
    })
  })
  
  cat("✅ PREFLIGHT PASSED - Safe to run full pipeline\n")
  invisible(TRUE)
}

#' Extract SCENIC metadata for enhanced edge annotation
#' 
#' @param scenicOptions SCENIC options object
#' @param scenicEdges data.frame of TF-target edges from SCENIC
#' @param verbose logical, whether to print status messages
#' @return data.frame with enhanced edge metadata
extractScenicMetadata <- function(scenicOptions, scenicEdges, verbose = FALSE) {
  
  if (verbose) {
    message("Extracting SCENIC metadata for enhanced edge annotation...")
  }
  
  # 1) Load artifacts ----------------------------------------------------------
  allTFs <- getDbTfs(scenicOptions)
  if (length(allTFs) < 100) {
    stop("CRITICAL FAILURE: Too few TFs found in SCENIC database (", length(allTFs), "). Check SCENIC setup.")
  }
  
  # Load regulonTargetsInfo as authoritative source of TF->target relationships
  if (verbose) message("Loading regulonTargetsInfo (authoritative TF->target source)...")
  tryCatch({
    regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
  }, error = function(e) {
    stop("CRITICAL FAILURE: Cannot load regulonTargetsInfo: ", e$message)
  })
  
  # Load motifEnrichment for annotation only (not for creating edges)
  if (verbose) {
    me_path <- tryCatch(getIntName(scenicOptions, "motifEnrichment"),
                        error = function(e) "")
    message("Loading motifEnrichment (for annotation only)... (", me_path, ")")
  }
  tryCatch({
    motifEnrichment <- loadInt(scenicOptions, "motifEnrichment")
  }, error = function(e) {
    if (verbose) message("Motif enrichment file not configured; proceeding with RI-provided motif evidence only.")
    motifEnrichment <- NULL
  })
  
  # Extract GENIE3 adjacencies
  genie3Links <- NULL
  genie3Keys <- c("genie3ll", "genie3wm", "genie3links", "adjacencies")
  
  for (key in genie3Keys) {
    fileName <- getIntName(scenicOptions, key)
    if (verbose) message("Attempting to load GENIE3 key '", key, "' from file: ", fileName)
    
    tryCatch({
      genie3Links <- loadInt(scenicOptions, key)
      if (!is.null(genie3Links)) {
        if (verbose) message("Found GENIE3 data with key: ", key)
        
        # Handle different GENIE3 formats
        if (key == "genie3wm" && is.matrix(genie3Links)) {
          if (verbose) message("Converting GENIE3 weight matrix to link list format")
          genie3Links <- as.data.frame(as.table(genie3Links))
          colnames(genie3Links) <- c("from", "to", "importance")
          genie3Links <- genie3Links[genie3Links$importance > 0, ]
        }
        break
      }
    }, error = function(e) {
      if (verbose) message("Failed to load key '", key, "': ", e$message)
    })
  }
  
  if (is.null(genie3Links)) {
    warning("No GENIE3 adjacency data found. Tried keys: ", paste(genie3Keys, collapse = ", "))
    genie3Links <- data.frame(from = character(0), to = character(0), importance = numeric(0))
  }
  
  # 2) Validate columns --------------------------------------------------------
  needColsRI <- c("TF", "target")
  
  # Handle SCENIC column naming variation: some versions use 'gene' instead of 'target'
  if ("gene" %in% colnames(regulonTargetsInfo) && "target" %in% colnames(regulonTargetsInfo)) {
    if (verbose) message("Both 'gene' and 'target' present in regulonTargetsInfo; using 'target' and dropping 'gene'.")
    regulonTargetsInfo$gene <- NULL
  } else if ("gene" %in% colnames(regulonTargetsInfo) && !"target" %in% colnames(regulonTargetsInfo)) {
    if (verbose) message("Converting regulonTargetsInfo 'gene' column to 'target'")
    names(regulonTargetsInfo)[names(regulonTargetsInfo) == "gene"] <- "target"
  }
  
  if (!all(needColsRI %in% colnames(regulonTargetsInfo))) {
    stop("CRITICAL FAILURE: regulonTargetsInfo missing required columns: ", 
         paste(setdiff(needColsRI, colnames(regulonTargetsInfo)), collapse = ", "),
         ". Found columns: ", paste(colnames(regulonTargetsInfo), collapse = ", "))
  }
  
  if (!is.null(motifEnrichment)) {
    needColsME <- c("motif", "geneSet")
    if (!all(needColsME %in% colnames(motifEnrichment))) {
      warning("motifEnrichment missing expected columns; motif annotation will be limited.")
      motifEnrichment <- NULL
    }
  }
  
  # 3) Canonical edges from regulonTargetsInfo only ----------------------------
  edgesBase <- regulonTargetsInfo %>%
    transmute(
      TF = .data$TF,
      Target = .data$target
    ) %>%
    distinct()
  
  if (nrow(edgesBase) == 0) {
    stop("CRITICAL FAILURE: No TF->target edges found in regulonTargetsInfo.")
  }
  
  if (verbose) {
    message("Found ", nrow(edgesBase), " canonical TF->target edges from regulonTargetsInfo")
  }
  
  # 4) Attach GENIE3 importance and validate correlation -----------------------
  if (!all(c("TF", "Target") %in% colnames(scenicEdges))) {
    stop("CRITICAL FAILURE: scenicEdges must have TF, Target columns before annotation.")
  }
  
  # Ensure correlation column exists with safe defaults
  if (!"corr" %in% colnames(scenicEdges)) {
    if (verbose) message("No 'corr' column found in scenicEdges, using default value 0")
    scenicEdges$corr <- 0
  }
  # Clean any non-finite correlation values
  scenicEdges$corr[!is.finite(scenicEdges$corr)] <- 0
  
  # Normalize GENIE3 data column names
  if (nrow(genie3Links) > 0) {
    if ("from" %in% colnames(genie3Links) && "to" %in% colnames(genie3Links)) {
      if ("importance" %in% colnames(genie3Links)) {
        genie3Links <- genie3Links %>% dplyr::rename(TF = from, Target = to, weight = importance)
      } else {
        stop("CRITICAL FAILURE: GENIE3 data has from/to columns but missing importance/weight.")
      }
    }
    
    expectedCols <- c("TF", "Target", "weight")
    missingCols <- setdiff(expectedCols, colnames(genie3Links))
    if (length(missingCols) > 0) {
      stop("CRITICAL FAILURE: GENIE3 data missing required columns: ", paste(missingCols, collapse=", "))
    }
    
    if (!is.numeric(genie3Links$weight)) {
      stop("CRITICAL FAILURE: GENIE3 weight column is not numeric.")
    }
    
    if (all(genie3Links$weight == 0, na.rm = TRUE)) {
      stop("CRITICAL FAILURE: All GENIE3 weights are zero - indicates failed GENIE3 inference.")
    }
  }
  
  # Join with GENIE3 data
  edges <- scenicEdges %>%
    right_join(edgesBase, by = c("TF", "Target")) %>%
    left_join(genie3Links %>% select(TF, Target, weight), by = c("TF", "Target")) %>%
    mutate(
      genie3Importance = ifelse(is.na(weight), 0, weight)
    )
  
  # 5) Optional: per-edge motif support from regulonTargetsInfo ---------------
  motifColsInRI <- intersect(colnames(regulonTargetsInfo), 
                             c("motif", "motifEnrichment", "motifNES", "motifAnnot", "NES", "AUC"))
  if (length(motifColsInRI) > 0) {
    if (verbose) message("Adding motif support from regulonTargetsInfo columns: ", paste(motifColsInRI, collapse=", "))
    motifSupport <- regulonTargetsInfo %>%
      dplyr::select(TF, target, all_of(motifColsInRI)) %>%
      dplyr::rename(Target = target) %>%
      dplyr::distinct()
    
    # Belt-and-suspenders: ensure join doesn't duplicate rows
    beforeJoin <- nrow(edges)
    edges <- edges %>% 
      left_join(motifSupport, by = c("TF", "Target"))
    stopifnot(nrow(edges) == beforeJoin)  # Guard against accidental duplication
    
    # Create motifNES if available
    if ("NES" %in% colnames(edges)) {
      edges <- edges %>% mutate(motifNES = ifelse(is.na(NES), 0, NES))
    }
  } else {
    # No per-edge motif data available from regulonTargetsInfo
    # Note: 0 here means 'no motif evidence provided', not evidence of absence
    if (verbose) message("No per-edge motif columns found in regulonTargetsInfo")
  }
  
  # ---- Ensure motif columns exist and are numeric-safe ----
  if (!"motifNES" %in% colnames(edges)) {
    edges$motifNES <- 0  # Recycled to length nrow(edges)
  }
  # If it exists but has NAs or non-finite values, sanitize:
  if (!is.numeric(edges$motifNES)) {
    # Try to coerce; if it fails, zero it out
    suppressWarnings({
      edges$motifNES <- as.numeric(edges$motifNES)
    })
  }
  edges$motifNES[!is.finite(edges$motifNES)] <- 0
  edges$motifNES[is.na(edges$motifNES)] <- 0
  
  # Create derived motif scores
  edges <- edges %>%
    mutate(
      motifConfidence = pmax(0, pmin(1, (motifNES - 2.0) / 3.0)),
      nMotifs = ifelse(motifNES > 0, 1, 0),  # Simple binary indicator
      
      # Compute overall prior strength with safe denominators
      priorStrength = ({
        weight95 <- quantile(genie3Importance, 0.95, na.rm = TRUE)
        if (is.na(weight95) || is.infinite(weight95) || weight95 <= 0) {
          weight95 <- max(genie3Importance, na.rm = TRUE)
          if (is.na(weight95) || weight95 <= 0) {
            weight95 <- 1e-9
          }
        }
        
        0.4 * pmax(0, pmin(1, genie3Importance / pmax(weight95, 1e-9))) +
          0.4 * motifConfidence +
          0.2 * pmax(0, pmin(1, nMotifs / 5))
      })
    )
  
  # 6) Final integrity checks --------------------------------------------------
  # Check TF recognition
  unrecognizedTFs <- setdiff(unique(edges$TF), allTFs)
  if (length(unrecognizedTFs) > 0) {
    recognitionRate <- 1 - (length(unrecognizedTFs) / length(unique(edges$TF)))
    if (recognitionRate < 0.8) {
      stop("CRITICAL FAILURE: Low TF recognition rate (", round(recognitionRate*100,1), 
           "%). Unrecognized TFs: ", paste(head(unrecognizedTFs, 5), collapse=", "))
    }
  }
  
  # Check required columns
  requiredCols <- c("TF", "Target", "genie3Importance", "motifConfidence", "priorStrength")
  missingCols <- setdiff(requiredCols, colnames(edges))
  if (length(missingCols) > 0) {
    stop("CRITICAL FAILURE: Enhanced edges missing required columns: ", paste(missingCols, collapse=", "))
  }
  
  # Check numeric validity
  if (any(!is.finite(edges$genie3Importance))) {
    stop("CRITICAL FAILURE: Non-finite genie3Importance detected.")
  }
  
  if (any(!is.finite(edges$priorStrength))) {
    stop("CRITICAL FAILURE: Non-finite priorStrength detected.")
  }
  
  if (verbose) {
    message("Successfully enhanced ", nrow(edges), " edges with SCENIC metadata")
    message("Validation passed: all edges have valid TF/Target pairs and numeric scores")
  }
  
  return(edges)
}

#' Assign edge signs using motif-aware approach
#' 
#' @param edges data.frame with TF, Target, corr, and SCENIC metadata
#' @param config configuration object with thresholds
#' @param verbose logical, whether to print status messages
#' @return data.frame with regType column assigned
assignMotifAwareSigns <- function(edges, config, verbose = FALSE) {
  
  if (verbose) {
    message("Assigning edge signs using motif-aware approach...")
  }
  
  # Motif-aware signing logic using config parameters
  edges <- edges %>%
    mutate(
      regType = case_when(
        # High-confidence motif predictions: trust them if correlation aligns
        motifConfidence > config$grnHighMotifThreshold & corr > config$grnMotifCorrAlignThreshold ~ "Activation",
        motifConfidence > config$grnHighMotifThreshold & corr < -config$grnMotifCorrAlignThreshold ~ "Inhibition",
        
        # Medium confidence: bias toward activation if positive correlation
        motifConfidence > config$grnMediumMotifThreshold & corr > config$grnPositiveThreshold ~ "Activation", 
        motifConfidence > config$grnMediumMotifThreshold & corr < -config$grnNegativeThreshold ~ "Inhibition",
        
        # Low motif support: fall back to correlation with config thresholds
        corr > config$grnPositiveThreshold ~ "Activation",
        corr < -config$grnNegativeThreshold ~ "Inhibition",
        
        # Default case
        TRUE ~ ifelse(corr >= 0, "Activation", "Inhibition")
      ),
      
      # Track sign confidence with more detail
      signConfidence = case_when(
        motifConfidence > config$grnHighMotifThreshold ~ "High_Motif",
        motifConfidence > config$grnMediumMotifThreshold ~ "Medium_Motif", 
        abs(corr) > config$grnPositiveThreshold ~ "High_Corr",
        TRUE ~ "Low_Confidence"
      ),
      
      # Compute numeric sign score for tie-breaking and diagnostics
      signScore = case_when(
        regType == "Activation" ~ pmax(motifConfidence, abs(corr)),
        regType == "Inhibition" ~ -pmax(motifConfidence, abs(corr)),
        TRUE ~ 0
      )
    )
  
  if (verbose) {
    # Handle potential list columns before counting
    tryCatch({
      signSummary <- edges %>% 
        # Ensure columns are vectors, not lists
        mutate(
          regType = as.character(regType),
          signConfidence = as.character(signConfidence)
        ) %>%
        count(regType, signConfidence) %>%
        arrange(regType, signConfidence)
      
      message("Sign assignment summary:")
      print(signSummary)
    }, error = function(e) {
      message("Could not generate sign summary: ", e$message)
      message("regType values: ", paste(unique(edges$regType), collapse = ", "))
      message("signConfidence values: ", paste(unique(edges$signConfidence), collapse = ", "))
    })
  }
  
  return(edges)
}

#' Compute composite regulator ranking scores
#' 
#' @param edges data.frame with correlation and SCENIC metadata
#' @param config configuration object with weights
#' @param verbose logical, whether to print status messages
#' @return data.frame with compositeScore column
computeCompositeRanking <- function(edges, config, verbose = FALSE) {
  
  if (verbose) {
    message("Computing composite regulator ranking scores...")
  }
  
  # Validate that weights sum to 1.0
  totalWeight <- config$grnCorrWeight + config$grnGenie3Weight + config$grnMotifWeight
  if (abs(totalWeight - 1.0) > 1e-6) {
    warning("Composite ranking weights do not sum to 1.0 (sum = ", round(totalWeight, 4), 
            "). Results may be uninterpretable.")
  }
  
  # Normalize scores to [0,1] range
  edges <- edges %>%
    mutate(
      # Normalize correlation by absolute value with epsilon protection
      corrScore = abs(corr) / pmax(max(abs(corr), na.rm = TRUE), 1e-9),
      
      # Normalize GENIE3 importance with epsilon protection
      genie3Score = genie3Importance / pmax(max(genie3Importance, na.rm = TRUE), 1e-9),
      
      # motifConfidence is already [0,1]
      
      # Compute composite score
      compositeScore = (
        config$grnCorrWeight * corrScore +
          config$grnGenie3Weight * genie3Score + 
          config$grnMotifWeight * motifConfidence
      )
    ) %>%
    # Compute rank within each target gene (grouped ranking only)
    group_by(Target) %>%
    mutate(regulatorRank = dplyr::dense_rank(dplyr::desc(compositeScore))) %>%
    ungroup()
  
  if (verbose) {
    message("Composite scoring completed. Top regulators by composite score:")
    topRegs <- edges %>%
      arrange(dplyr::desc(compositeScore)) %>%
      select(TF, Target, compositeScore, corrScore, genie3Score, motifConfidence) %>%
      head(10)
    print(topRegs)
  }
  
  return(edges)
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
createAndSaveGrnPlot <- function(graph, plotPath, title = "Gene Regulatory Network", 
                                 config, edges = NULL, verbose = FALSE) {
  
  # Validate inputs
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (igraph::vcount(graph) == 0) {
    warning("Graph has no vertices - creating empty plot")
  }
  
  # Create plot
  if (verbose) {
    message("Creating enhanced GRN plot with ", igraph::vcount(graph), " nodes and ", igraph::ecount(graph), " edges")
  }
  
  # Add enhanced edge information if available by pushing to graph attributes
  if (!is.null(edges)) {
    edgeList <- igraph::as_data_frame(graph, what = "edges")
    enhancedEdges <- edgeList %>%
      left_join(edges %>% select(from = TF, to = Target, regType, priorStrength, signConfidence), 
                by = c("from", "to"))
    
    # Push attributes to graph edges
    if ("regType" %in% colnames(enhancedEdges)) {
      igraph::E(graph)$regType <- enhancedEdges$regType
    }
    if ("priorStrength" %in% colnames(enhancedEdges)) {
      igraph::E(graph)$priorStrength <- enhancedEdges$priorStrength
    }
    if ("signConfidence" %in% colnames(enhancedEdges)) {
      igraph::E(graph)$signConfidence <- enhancedEdges$signConfidence
    }
  }
  
  # Create the plot (now attributes are accessible via graph)
  grnPlot <- ggraph::ggraph(graph, layout = config$networkLayout) +
    {
      if (igraph::ecount(graph) > 0) {
        if ("regType" %in% igraph::edge_attr_names(graph)) {
          if ("priorStrength" %in% igraph::edge_attr_names(graph)) {
            # Enhanced edges with prior strength
            ggraph::geom_edge_link(aes(colour = regType, alpha = priorStrength), 
                                   arrow = grid::arrow(length = grid::unit(3, "mm")))
          } else {
            # Standard edges with regType
            ggraph::geom_edge_link(aes(colour = regType), alpha = config$plotAlpha,
                                   arrow = grid::arrow(length = grid::unit(3, "mm")))
          }
        } else {
          # Basic edges
          ggraph::geom_edge_link(alpha = config$plotAlpha,
                                 arrow = grid::arrow(length = grid::unit(3, "mm")))
        }
      } else {
        ggraph::geom_edge_link(alpha = config$plotAlpha)
      }
    } +
    ggraph::geom_node_point(size = config$pointSize, color = "steelblue") +
    ggraph::geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    ggraph::scale_edge_colour_manual(values = c("Activation" = "darkgreen", 
                                                "Inhibition" = "darkred")) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = config$hjust)) +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_blank()
    )
  
  # Save plot
  if (verbose) {
    message("Saving enhanced GRN plot to: ", plotPath)
  }
  
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
createAndSaveGraphml <- function(graph, graphmlPath, edges = NULL, verbose = FALSE) {
  
  # Validate inputs
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (igraph::vcount(graph) == 0) {
    stop("Cannot save GraphML: graph has no vertices")
  }
  
  if (igraph::ecount(graph) == 0) {
    stop("Cannot save GraphML: graph has no edges")
  }
  
  if (verbose) {
    message("Saving GraphML file to: ", graphmlPath)
  }
  
  # Add edge attributes if available (validation-only, no inference)
  if (!is.null(edges)) {
    # Validate required columns exist
    requiredCols <- c("TF", "Target", "genie3Importance")
    missingCols <- setdiff(requiredCols, colnames(edges))
    if (length(missingCols) > 0) {
      stop("CRITICAL FAILURE: edges missing required columns for GraphML: ", paste(missingCols, collapse=", "))
    }
    
    edgeList <- igraph::as_data_frame(graph, what = "edges")
    enhancedEdges <- edgeList %>%
      left_join(edges %>% 
                select(from = TF, to = Target, genie3Importance, 
                       any_of(c("regType", "priorStrength", "signConfidence", "signScore", "motifNES"))), 
                by = c("from", "to"))
    
    # Add validated attributes back to graph (replace NA with defaults)
    availableCols <- intersect(c("genie3Importance", "regType", "priorStrength", "signConfidence", "signScore", "motifNES"), 
                               colnames(enhancedEdges))
    
    for (col in availableCols) {
      if (col %in% colnames(enhancedEdges)) {
        if (is.numeric(enhancedEdges[[col]])) {
          # Replace NA with 0 for numeric columns
          enhancedEdges[[col]][is.na(enhancedEdges[[col]])] <- 0
        } else {
          # Replace NA with "Unknown" for character columns
          enhancedEdges[[col]][is.na(enhancedEdges[[col]])] <- "Unknown"
        }
        igraph::E(graph)[[col]] <- enhancedEdges[[col]]
      }
    }
  }
  
  # Ensure directory exists
  dir.create(dirname(graphmlPath), recursive = TRUE, showWarnings = FALSE)
  
  # Save GraphML file
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
#' @return list with paths to saved files
saveEnhancedGrnOutputSet <- function(graph, edges, rdsPath, plotPath, graphmlPath, 
                                     edgesPath, title, config, verbose = FALSE) {
  
  if (verbose) {
    message("Saving complete enhanced GRN output set...")
  }
  
  # Save RDS file
  if (verbose) {
    message("Saving GRN RDS to: ", rdsPath)
  }
  saveRDS(graph, file = rdsPath)
  
  # Save enhanced edges
  if (verbose) {
    message("Saving enhanced edges to: ", edgesPath)
  }
  write.csv(edges, file = edgesPath, row.names = FALSE)
  
  # Save enhanced plot
  createAndSaveGrnPlot(graph, plotPath, title, config, edges, verbose)
  
  # Save enhanced GraphML
  createAndSaveGraphml(graph, graphmlPath, edges, verbose)
  
  # Return paths for verification
  savedFiles <- list(
    rds = rdsPath,
    edges = edgesPath,
    plot = plotPath,
    graphml = graphmlPath
  )
  
  if (verbose) {
    message("Completed saving enhanced GRN output set")
  }
  
  invisible(savedFiles)
}