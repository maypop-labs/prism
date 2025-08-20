# =============================================================================
# grnManager.R (Enhanced & Cleaned)
# Functions for creating, visualizing, and saving gene regulatory networks
# with improved integration of SCENIC/GENIE3 biological priors
# =============================================================================

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
  
  # Load SCENIC intermediate results
  tryCatch({
    # Extract motif enrichment data
    motifEnrichment <- loadInt(scenicOptions, "motifEnrichment")
    
    # Extract regulon target information
    regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo") 
    
    # Extract GENIE3 adjacencies with fallback for different SCENIC versions
    genie3Links <- NULL
    genie3Keys <- c("genie3ll", "genie3links", "adjacencies")
    
    for (key in genie3Keys) {
      tryCatch({
        genie3Links <- loadInt(scenicOptions, key)
        if (!is.null(genie3Links)) {
          if (verbose) message("Found GENIE3 data with key: ", key)
          break
        }
      }, error = function(e) NULL)
    }
    
    if (is.null(genie3Links)) {
      warning("No GENIE3 adjacency data found. Tried keys: ", paste(genie3Keys, collapse = ", "))
      genie3Links <- data.frame(from = character(0), to = character(0), importance = numeric(0))
    }
    
    if (verbose) {
      message("Found SCENIC metadata: motif enrichment (", nrow(motifEnrichment), 
              " entries), regulon targets (", nrow(regulonTargetsInfo), " entries)")
    }
    
    # Prepare motif data for merging
    motifData <- motifEnrichment %>%
      group_by(TF, target) %>%
      summarize(
        motifNES = max(NES, na.rm = TRUE),
        motifAUC = max(AUC, na.rm = TRUE), 
        nMotifs = n(),
        .groups = "drop"
      ) %>%
      rename(Target = target)
    
    # Prepare GENIE3 importance scores
    genie3Data <- genie3Links %>%
      rename(TF = from, Target = to, genie3Importance = importance)
    
    # Merge with edge data
    enhancedEdges <- scenicEdges %>%
      left_join(motifData, by = c("TF", "Target")) %>%
      left_join(genie3Data, by = c("TF", "Target")) %>%
      mutate(
        # Fill missing values with defaults (explicit handling)
        motifNES = ifelse(is.na(motifNES) | is.infinite(motifNES), 0, motifNES),
        motifAUC = ifelse(is.na(motifAUC) | is.infinite(motifAUC), 0, motifAUC),
        nMotifs = ifelse(is.na(nMotifs), 0, nMotifs),
        genie3Importance = ifelse(is.na(genie3Importance), 0, genie3Importance),
        
        # Compute motif confidence score (explicit 0 for missing)
        motifConfidence = ifelse(is.na(motifNES), 0, pmax(0, pmin(1, (motifNES - 2.0) / 3.0))),
        
        # Compute overall prior strength with safe denominators
        priorStrength = (
          0.4 * pmax(0, pmin(1, genie3Importance / pmax(quantile(genie3Importance, 0.95, na.rm = TRUE), 1e-9))) +
            0.4 * motifConfidence +
            0.2 * pmax(0, pmin(1, nMotifs / 5))
        )
      )
    
    if (verbose) {
      message("Enhanced ", nrow(enhancedEdges), " edges with SCENIC metadata")
    }
    
    return(enhancedEdges)
    
  }, error = function(e) {
    warning("Failed to extract SCENIC metadata: ", e$message, 
            ". Proceeding with basic edge data.")
    
    # Return original edges with default columns
    return(scenicEdges %>%
             mutate(
               motifNES = 0,
               motifAUC = 0,
               nMotifs = 0,
               genie3Importance = 0,
               motifConfidence = 0,
               priorStrength = 0
             ))
  })
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
    signSummary <- edges %>% 
      count(regType, signConfidence) %>%
      arrange(regType, signConfidence)
    
    message("Sign assignment summary:")
    print(signSummary)
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
      # Normalize correlation by absolute value
      corrScore = abs(corr) / max(abs(corr), na.rm = TRUE),
      
      # Normalize GENIE3 importance 
      genie3Score = genie3Importance / max(genie3Importance, na.rm = TRUE),
      
      # motifConfidence is already [0,1]
      
      # Compute composite score
      compositeScore = (
        config$grnCorrWeight * corrScore +
          config$grnGenie3Weight * genie3Score + 
          config$grnMotifWeight * motifConfidence
      ),
      
      # Compute rank within each target gene (deterministic)
      regulatorRank = dense_rank(desc(compositeScore))
    ) %>%
    group_by(Target) %>%
    mutate(regulatorRank = dense_rank(desc(compositeScore))) %>%
    ungroup()
  
  if (verbose) {
    message("Composite scoring completed. Top regulators by composite score:")
    topRegs <- edges %>%
      arrange(desc(compositeScore)) %>%
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
  
  if (verbose) {
    message("Saving enhanced GraphML file to: ", graphmlPath)
  }
  
  # Add edge attributes if available
  if (!is.null(edges)) {
    edgeList <- igraph::as_data_frame(graph, what = "edges")
    enhancedEdges <- edgeList %>%
      left_join(edges %>% select(from = TF, to = Target, regType, priorStrength, 
                                 signConfidence, signScore, motifNES, genie3Importance), 
                by = c("from", "to"))
    
    # Add attributes back to graph
    for (col in c("regType", "priorStrength", "signConfidence", "signScore", "motifNES", "genie3Importance")) {
      if (col %in% colnames(enhancedEdges)) {
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
      message("Successfully saved enhanced GraphML with ", igraph::vcount(graph), 
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