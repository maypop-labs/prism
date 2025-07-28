# =============================================================================
# setupManager.R
# Purpose: Configuration loading, library management, and script initialization
# Dependencies: yaml, doParallel, foreach
# =============================================================================

# =============================================================================
# Path Management
# =============================================================================

#' Build file paths for PRISM project outputs
#'
#' Constructs standardized file paths based on the project configuration and
#' optional cell type and trajectory parameters.
#'
#' @param type Directory type ("plots", "rds", "tsv", "txt", "monocle3", "graphml")
#' @param filename Optional filename to append
#' @param cellType Optional cell type for path construction
#' @param trajectory Optional trajectory for path construction
#' @return Character string with complete file path
#' @export
#' @examples
#' buildPath("rds", "merged_seurat.rds")
#' buildPath("plots", "figure1.png", cellType = "keratinocytes")
buildPath <- function(type, filename = NULL, cellType = NULL, trajectory = NULL) {
  
  # Validate type parameter
  validTypes <- c("plots", "rds", "tsv", "txt", "monocle3", "graphml")
  if (!type %in% validTypes) {
    stop("Invalid type: ", type, ". Must be one of: ", paste(validTypes, collapse = ", "))
  }
  
  # Get config (assumes global config variable exists)
  if (!exists("config", envir = .GlobalEnv)) {
    stop("Configuration not loaded. Call initializeScript() first.")
  }
  config <- get("config", envir = .GlobalEnv)
  
  # Build base path
  basePath <- paste0(config$rootPath, "results/", type, "/")
  
  # Return base path if no filename specified
  if (is.null(filename)) {
    return(basePath)
  }
  
  # Construct filename with optional prefixes
  if (!is.null(cellType) && !is.null(trajectory)) {
    fullFilename <- paste0(cellType, "_", trajectory, "_", filename)
  } else if (!is.null(cellType)) {
    fullFilename <- paste0(cellType, "_", filename)
  } else {
    fullFilename <- filename
  }
  
  return(paste0(basePath, fullFilename))
}

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
    "ggplot2", "ggraph", "data.table", "mgcv", "gmp", "progress", "progressr",
    
    # Extended analysis
    "tidyverse"
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
  requiredFields <- c("rootPath", "ageVec", "nameVec", "cores", "saveResults")
  missingFields <- setdiff(requiredFields, names(config))
  
  if (length(missingFields) > 0) {
    stop("Missing required configuration fields: ", paste(missingFields, collapse = ", "))
  }
  
  # Validate data consistency
  if (length(config$ageVec) != length(config$nameVec)) {
    stop("ageVec and nameVec must have the same length")
  }
  
  # Normalize root path
  if (!endsWith(config$rootPath, "/") && !endsWith(config$rootPath, "\\")) {
    config$rootPath <- paste0(config$rootPath, "/")
  }
  
  # Validate numeric parameters
  numericParams <- c("cores", "ageCorrelation", "boolMaxRegulators", 
                     "minCellTypeNumber", "figWidth", "figHeight", "figDPI")
  
  for (param in numericParams) {
    if (param %in% names(config) && !is.numeric(config[[param]])) {
      warning("Parameter '", param, "' should be numeric, got: ", class(config[[param]]))
    }
  }
  
  # Set reasonable defaults for missing optional parameters
  defaultValues <- list(
    cores = 1,
    ageCorrelation = 0.6,
    boolMaxRegulators = 3,
    minCellTypeNumber = 1000,
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
  
  # Configure Seurat version if specified
  if ("SeuratAssay" %in% names(config) && !is.null(config$SeuratAssay)) {
    options(Seurat.object.assay.version = config$SeuratAssay)
    message("Set Seurat assay version to: ", config$SeuratAssay)
  }
  
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
# Utility Functions
# =============================================================================

#' Time execution of R expressions with formatted output
#'
#' Wrapper function that times the execution of an R expression and reports
#' the elapsed time in a human-readable format.
#'
#' @param expr R expression to time
#' @param name Descriptive name for the operation (default: "task")
#' @return Result of the expression (timing is printed as side effect)
#' @export
#' @examples
#' result <- timeIt({
#'   Sys.sleep(2)
#'   1 + 1
#' }, "test calculation")
timeIt <- function(expr, name = "task") {
  
  startTime <- Sys.time()
  result <- force(expr)
  endTime <- Sys.time()
  
  # Calculate duration
  duration <- difftime(endTime, startTime, units = "secs")
  durationSeconds <- as.numeric(duration)
  
  # Format output
  if (durationSeconds >= 60) {
    minutes <- floor(durationSeconds / 60)
    seconds <- durationSeconds %% 60
    timeStr <- sprintf("%d min %.1f sec", minutes, seconds)
  } else {
    timeStr <- sprintf("%.1f sec", durationSeconds)
  }
  
  message(sprintf("%s completed in %s", name, timeStr))
  
  invisible(result)
}

#' Check system requirements and provide diagnostic information
#'
#' Examines the current system environment and reports information relevant
#' to PRISM pipeline execution, including memory, disk space, and R version.
#'
#' @param config Configuration list from loadProjectConfig()
#' @export
checkSystemRequirements <- function(config) {
  
  message("=== System Requirements Check ===")
  
  # Check R version
  rVersion <- R.version.string
  message("R version: ", rVersion)
  
  if (getRversion() < "4.0.0") {
    warning("R version < 4.0 may cause compatibility issues with some packages")
  }
  
  # Check available memory (Unix systems only)
  if (.Platform$OS.type == "unix") {
    memInfo <- try(system("free -h", intern = TRUE, ignore.stderr = TRUE), silent = TRUE)
    if (!inherits(memInfo, "try-error") && length(memInfo) >= 2) {
      message("System memory information:")
      cat(memInfo[1:2], sep = "\n")
    }
  } else {
    message("Memory check not available on this platform")
  }
  
  # Check output directory
  if (dir.exists(config$rootPath)) {
    message("Output directory exists: ", config$rootPath)
  } else {
    message("Output directory will be created: ", config$rootPath)
  }
  
  # Check disk space (if possible)
  tryCatch({
    diskInfo <- file.info(dirname(config$rootPath))
    if (!is.na(diskInfo$size)) {
      message("Output directory parent accessible")
    }
  }, error = function(e) {
    warning("Could not access output directory: ", e$message)
  })
  
  # Check critical package availability
  criticalPackages <- c("Seurat", "monocle3", "BoolNet", "yaml")
  missingCritical <- character()
  
  for (pkg in criticalPackages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missingCritical <- c(missingCritical, pkg)
    }
  }
  
  if (length(missingCritical) > 0) {
    warning("Critical packages not installed: ", paste(missingCritical, collapse = ", "))
    message("Install with: install.packages(c('", paste(missingCritical, collapse = "', '"), "'))")
  } else {
    message("All critical packages are available")
  }
  
  # Check parallel processing capability
  if (config$cores > 1) {
    availableCores <- parallel::detectCores()
    if (is.na(availableCores)) {
      warning("Cannot detect number of available cores")
    } else if (config$cores > availableCores) {
      warning("Requested cores (", config$cores, ") exceeds available cores (", availableCores, ")")
    } else {
      message("Parallel processing configured: ", config$cores, " of ", availableCores, " cores")
    }
  }
  
  message("System check complete\n")
  invisible(TRUE)
}

# =============================================================================
# File Validation
# =============================================================================

#' Validate that required files exist before script execution
#'
#' Checks for the existence of required input files and provides clear error
#' messages if files are missing, helping users understand dependencies.
#'
#' @param ... File paths to validate (can be named for better error messages)
#' @param stopOnMissing Whether to stop execution if files are missing (default: TRUE)
#' @return Named logical vector indicating which files exist
#' @export
#' @examples
#' validateRequiredFiles(
#'   "merged seurat" = "results/rds/merged_seurat.rds",
#'   "cell types" = "results/rds/cell_types.rds"
#' )
validateRequiredFiles <- function(..., stopOnMissing = TRUE) {
  
  files <- list(...)
  
  if (length(files) == 0) {
    return(logical(0))
  }
  
  # Handle both named and unnamed arguments
  fileNames <- names(files)
  if (is.null(fileNames)) {
    fileNames <- as.character(files)
  } else {
    # Replace empty names with file paths
    emptyNames <- fileNames == ""
    fileNames[emptyNames] <- as.character(files)[emptyNames]
  }
  
  # Check file existence
  exists <- sapply(files, file.exists)
  names(exists) <- fileNames
  
  missing <- !exists
  
  if (any(missing)) {
    missingFiles <- names(exists)[missing]
    errorMsg <- paste0(
      "Required files not found:\n",
      paste(paste("  -", missingFiles), collapse = "\n")
    )
    
    if (stopOnMissing) {
      stop(errorMsg, "\n\nEnsure prerequisite scripts have been run successfully.")
    } else {
      warning(errorMsg)
    }
  } else {
    message("All required files validated successfully")
  }
  
  return(exists)
}

# =============================================================================
# Progress Logging
# =============================================================================

#' Clear console for clean output display
#'
#' Clears the R console if running interactively, providing clean output
#' display for script execution.
#'
#' @export
clearConsole <- function() {
  if (interactive()) {
    cat("\014\n")
  }
}

#' Log progress messages with optional step counting
#'
#' Provides standardized progress logging with optional step numbering
#' for clear tracking of script execution.
#'
#' @param message Message to display
#' @param step Current step number (optional)
#' @param total Total number of steps (optional)
#' @export
#' @examples
#' logProgress("Loading data")
#' logProgress("Processing genes", step = 5, total = 100)
logProgress <- function(message, step = NULL, total = NULL) {
  if (!is.null(step) && !is.null(total)) {
    cat(sprintf("[%d/%d] %s\n", step, total, message))
  } else {
    message(message)
  }
}

#' Log section headers for script organization
#'
#' Creates visually distinct section headers to organize script output
#' and improve readability.
#'
#' @param title Section title
#' @export
#' @examples
#' logSection("Data Loading and QC")
logSection <- function(title) {
  message("=== ", title, " ===")
}