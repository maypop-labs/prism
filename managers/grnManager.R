# =============================================================================
# grnManager.R
# Functions for creating, visualizing, and saving gene regulatory networks
# =============================================================================

#' Create and save a GRN plot
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param plotPath character, full path where plot should be saved
#' @param title character, plot title (default: "Gene Regulatory Network")
#' @param config list, configuration object containing plot parameters
#' @param verbose logical, whether to print status messages
#' @return ggplot object (invisibly)
createAndSaveGrnPlot <- function(graph, plotPath, title = "Gene Regulatory Network", config, verbose = FALSE) {
  
  # Validate inputs
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (igraph::vcount(graph) == 0) {
    warning("Graph has no vertices - creating empty plot")
  }
  
  # Create plot
  if (verbose) {
    message("Creating GRN plot with ", igraph::vcount(graph), " nodes and ", igraph::ecount(graph), " edges")
  }
  
  # Handle edge attributes safely
  edgeData <- igraph::as_data_frame(graph, what = "edges")
  
  # Create the plot
  grnPlot <- ggraph::ggraph(graph, layout = config$networkLayout) +
    {
      if (nrow(edgeData) > 0 && "regType" %in% colnames(edgeData)) {
        ggraph::geom_edge_link(aes(colour = regType), alpha = config$plotAlpha)
      } else {
        ggraph::geom_edge_link(alpha = config$plotAlpha)
      }
    } +
    ggraph::geom_node_point(size = config$pointSize) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = config$hjust)) +
    ggraph::theme_graph(background = "white") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
  # Save plot
  if (verbose) {
    message("Saving GRN plot to: ", plotPath)
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

#' Create and save a GraphML file
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param graphmlPath character, full path where GraphML file should be saved
#' @param verbose logical, whether to print status messages
#' @return NULL (invisibly)
createAndSaveGraphml <- function(graph, graphmlPath, verbose = FALSE) {
  
  # Validate inputs
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  if (verbose) {
    message("Saving GraphML file to: ", graphmlPath)
  }
  
  # Ensure directory exists
  dir.create(dirname(graphmlPath), recursive = TRUE, showWarnings = FALSE)
  
  # Save GraphML file
  tryCatch({
    igraph::write_graph(graph, graphmlPath, format = "graphml")
    if (verbose) {
      message("Successfully saved GraphML with ", igraph::vcount(graph), " nodes and ", igraph::ecount(graph), " edges")
    }
  }, error = function(e) {
    stop("Failed to save GraphML file: ", e$message)
  })
  
  invisible(NULL)
}

#' Save complete GRN output set (RDS + plot + GraphML)
#' 
#' @param graph igraph object representing the gene regulatory network
#' @param rdsPath character, path for RDS file
#' @param plotPath character, path for plot file
#' @param graphmlPath character, path for GraphML file
#' @param title character, plot title
#' @param config list, configuration object
#' @param verbose logical, whether to print status messages
#' @return list with paths to saved files
saveGrnOutputSet <- function(graph, rdsPath, plotPath, graphmlPath, title, config, verbose = FALSE) {
  
  if (verbose) {
    message("Saving complete GRN output set...")
  }
  
  # Save RDS file
  if (verbose) {
    message("Saving GRN RDS to: ", rdsPath)
  }
  saveRDS(graph, file = rdsPath)
  
  # Save plot
  createAndSaveGrnPlot(graph, plotPath, title, config, verbose)
  
  # Save GraphML
  createAndSaveGraphml(graph, graphmlPath, verbose)
  
  # Return paths for verification
  savedFiles <- list(
    rds = rdsPath,
    plot = plotPath,
    graphml = graphmlPath
  )
  
  if (verbose) {
    message("Completed saving GRN output set")
  }
  
  invisible(savedFiles)
}

#' Get network summary statistics
#' 
#' @param graph igraph object
#' @param label character, label for the network (e.g., "preprocessed", "final")
#' @return named list of network statistics
getNetworkStats <- function(graph, label = "") {
  
  if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be an igraph object")
  }
  
  stats <- list(
    label = label,
    nodes = igraph::vcount(graph),
    edges = igraph::ecount(graph),
    density = igraph::edge_density(graph),
    avgDegree = ifelse(igraph::vcount(graph) > 0, mean(igraph::degree(graph)), 0),
    isConnected = igraph::is.connected(graph, mode = "weak"),
    numComponents = igraph::components(graph)$no,
    numStrongComponents = igraph::components(graph, mode = "strong")$no
  )
  
  return(stats)
}

#' Print network comparison
#' 
#' @param graphBefore igraph object (before processing)
#' @param graphAfter igraph object (after processing)
#' @param verbose logical, whether to print detailed comparison
printNetworkComparison <- function(graphBefore, graphAfter, verbose = TRUE) {
  
  if (!verbose) return(invisible(NULL))
  
  statsBefore <- getNetworkStats(graphBefore, "before")
  statsAfter <- getNetworkStats(graphAfter, "after")
  
  message("\n=== Network Processing Summary ===")
  message("Nodes: ", statsBefore$nodes, " -> ", statsAfter$nodes, 
          " (", statsAfter$nodes - statsBefore$nodes, ")")
  message("Edges: ", statsBefore$edges, " -> ", statsAfter$edges, 
          " (", statsAfter$edges - statsBefore$edges, ")")
  message("Density: ", round(statsBefore$density, 4), " -> ", round(statsAfter$density, 4))
  message("Components: ", statsBefore$numComponents, " -> ", statsAfter$numComponents)
  message("===================================\n")
  
  invisible(list(before = statsBefore, after = statsAfter))
}