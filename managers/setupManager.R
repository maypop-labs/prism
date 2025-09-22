# =============================================================================
# setupManager.R
# Purpose: Configuration loading, library management, and script initialization
# Dependencies: yaml, doParallel, foreach
# =============================================================================

# =============================================================================
# Library Management
# =============================================================================

#' Load required libraries for PRISM analysis
#'
#' Loads all standard libraries required for PRISM pipeline execution.
#' Provides error handling and reports any failed library loads.
#'
#' @param suppressMessages Whether to suppress package startup messages (default: TRUE)
#' @return Invisible TRUE if successful, stops execution if critical libraries fail
#' @export
loadRequiredLibraries <- function(suppressMessages = TRUE) {
  
  # Define required libraries
  standardLibs <- c(
    # Core R packages
    "foreach", "doParallel", "yaml", "tools", "dplyr", "tidyr", "purrr", "readr",

    
    # Single-cell analysis
    "Seurat", "SeuratData", "SeuratDisk", "monocle3", "SingleCellExperiment",
    "celldex", "SingleR", "Matrix",
    
    # Network analysis
    "BoolNet", "igraph", "SCENIC", "GENIE3", "GeneSwitches",
    
    # Visualization and utilities
    "ggplot2", "ggraph", "data.table", "mgcv", "gmp", "patchwork", "progress",
    "progressr", "gridExtra", "ComplexHeatmap", "jsonlite",
    
    # Extended analysis
    "tidyverse", "lme4", "lmerTest"
  )
  
  # Track loading results
  loadedLibs <- character()
  failedLibs <- character()
  
  message("Loading required libraries...")
  
  # Load each library with error handling
  for (lib in standardLibs) {
    tryCatch({
      if (suppressMessages) {
        suppressPackageStartupMessages(library(lib, character.only = TRUE))
      } else {
        library(lib, character.only = TRUE)
      }
      loadedLibs <- c(loadedLibs, lib)
    }, error = function(e) {
      failedLibs <- c(failedLibs, lib)
      warning("Failed to load library '", lib, "': ", e$message)
    })
  }
  
  # Report results
  message("Successfully loaded ", length(loadedLibs), " libraries")
  
  if (length(failedLibs) > 0) {
    warning("Failed to load ", length(failedLibs), " libraries: ", 
            paste(failedLibs, collapse = ", "))
    message("\nTo install missing packages, run:")
    message("install.packages(c('", paste(failedLibs, collapse = "', '"), "'))")
    
    # Check for critical libraries
    criticalLibs <- c("Seurat", "monocle3", "BoolNet", "yaml", "dplyr")
    criticalMissing <- intersect(failedLibs, criticalLibs)
    
    if (length(criticalMissing) > 0) {
      stop("Critical libraries missing: ", paste(criticalMissing, collapse = ", "), 
           "\nInstall these packages before continuing.")
    }
  }
  
  invisible(TRUE)
}

# =============================================================================
# Configuration Management
# =============================================================================

#' Load and validate PRISM configuration
#'
#' Loads the project configuration from config.yaml and validates required
#' parameters for proper pipeline execution.
#'
#' @param configPath Path to configuration YAML file (default: "config.yaml")
#' @return Configuration list with validated parameters
#' @export
loadProjectConfig <- function(configPath = "config.yaml") {
  
  # Check if config file exists
  if (!file.exists(configPath)) {
    stop("Configuration file not found: ", configPath)
  }
  
  # Load configuration
  config <- yaml::read_yaml(configPath)
  
  # Validate required fields
  requiredFields <- c("rootPath", "ages", "donorIDs", "cores", "saveResults")
  missingFields <- setdiff(requiredFields, names(config))
  
  if (length(missingFields) > 0) {
    stop("Missing required configuration fields: ", paste(missingFields, collapse = ", "))
  }
  
  # Validate data consistency
  if (length(config$ages) != length(config$donorIDs)) {
    stop("ages and donorIDs must have the same length")
  }
  
  # Normalize root path
  if (!endsWith(config$rootPath, "/") && !endsWith(config$rootPath, "\\")) {
    config$rootPath <- paste0(config$rootPath, "/")
  }
  
  # Validate numeric parameters
  numericParams <- c("cores", "pseudotimeMinAgeCorrelation", "boolMaxRegulators", 
                     "singleRMinimumNumberOfCells", "figWidth", "figHeight", "figDPI")
  
  for (param in numericParams) {
    if (param %in% names(config) && !is.numeric(config[[param]])) {
      warning("Parameter '", param, "' should be numeric, got: ", class(config[[param]]))
    }
  }
  
  # Set reasonable defaults for missing optional parameters
  defaultValues <- list(
    cores = 1,
    pseudotimeMinAgeCorrelation = 0.6,
    boolMaxRegulators = 3,
    singleRMinimumNumberOfCells = 1000,
    saveResults = TRUE,
    verbose = TRUE,
    figWidth = 6.5,
    figHeight = 5,
    figDPI = 300
  )
  
  for (param in names(defaultValues)) {
    if (!param %in% names(config)) {
      config[[param]] <- defaultValues[[param]]
      message("Using default value for '", param, "': ", defaultValues[[param]])
    }
  }
  
  message("Configuration loaded successfully from: ", configPath)
  return(config)
}

# =============================================================================
# R Environment Setup
# =============================================================================

#' Configure R options for PRISM analysis
#'
#' Sets standard R options including warning suppression, Seurat configuration,
#' and parallel processing setup based on the loaded configuration.
#'
#' @param config Configuration list from loadProjectConfig()
#' @export
setupREnvironment <- function(config) {
  
  # Suppress warnings during processing (can be noisy with large datasets)
  options(warn = -1)
  
  # Configure parallel processing
  if ("cores" %in% names(config) && config$cores > 1) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      warning("doParallel package required for multi-core processing")
    } else {
      doParallel::registerDoParallel(cores = config$cores)
      message("Registered ", config$cores, " cores for parallel processing")
    }
  }
  
  # Set random seed for reproducibility if specified
  if ("randomSeed" %in% names(config)) {
    set.seed(config$randomSeed)
    message("Set random seed to: ", config$randomSeed)
  }
  
  # Configure memory management for large datasets
  if ("maxMemoryGB" %in% names(config)) {
    memoryBytes <- config$maxMemoryGB * 1024^3
    options(future.globals.maxSize = memoryBytes)
    message("Set maximum memory limit to: ", config$maxMemoryGB, " GB")
  }
  
  # Configure progress reporting
  if ("verbose" %in% names(config) && config$verbose) {
    options(progressr.enable = TRUE)
  }
  
  invisible(TRUE)
}

# =============================================================================
# Directory Management
# =============================================================================

#' Ensure all required directories exist
#'
#' Creates the standard directory structure for PRISM project outputs
#' if directories don't already exist.
#'
#' @param config Configuration list with rootPath
#' @export
ensureDirectories <- function(config) {
  
  # Define standard directory structure
  directories <- c(
    "results",
    "results/plots", 
    "results/rds",
    "results/tsv",
    "results/txt", 
    "results/monocle3",
    "results/graphml"
  )
  
  # Create each directory
  created <- character()
  for (dir in directories) {
    fullPath <- paste0(config$rootPath, dir)
    if (!dir.exists(fullPath)) {
      dir.create(fullPath, recursive = TRUE, showWarnings = FALSE)
      created <- c(created, dir)
    }
  }
  
  if (length(created) > 0) {
    message("Created directories: ", paste(created, collapse = ", "))
  }
  
  invisible(TRUE)
}

# =============================================================================
# Complete Script Initialization
# =============================================================================

#' Initialize PRISM script with complete setup
#'
#' Performs complete initialization including configuration loading, library
#' loading, R environment setup, and directory creation. This is the main
#' function that scripts should call for setup.
#'
#' @param configPath Path to configuration file (default: "config.yaml")
#' @return Configuration list for use in the calling script
#' @export
initializeScript <- function(configPath = "config.yaml") {
  
  message("=== Initializing PRISM Script ===")
  
  # Load configuration
  config <- loadProjectConfig(configPath)
  
  # Set up R environment
  setupREnvironment(config)
  
  # Load required libraries
  loadRequiredLibraries(suppressMessages = TRUE)
  
  # Create necessary directories
  ensureDirectories(config)
  
  # Store config globally for use by other functions
  assign("config", config, envir = .GlobalEnv)
  
  message("Script initialization complete")
  message("Ready to begin analysis\n")
  
  return(config)
}
