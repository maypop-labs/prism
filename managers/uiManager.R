# =============================================================================
# uiManager.R
# Purpose: User interface functions for interactive script elements
# Dependencies: None (base R only)
# =============================================================================

# =============================================================================
# Console Management
# =============================================================================

#' Clear console for clean output display
#'
#' Clears the R console if running interactively, providing a clean slate
#' for script output and user interaction.
#'
#' @export
#' @examples
#' clearConsole()  # Clears screen if running interactively
clearConsole <- function() {
  if (interactive()) {
    cat("\014\n")
  }
}

# =============================================================================
# Progress and Message Logging
# =============================================================================

#' Log progress messages with optional step counting
#'
#' Provides standardized progress logging with optional step numbering for
#' clear tracking of long-running operations.
#'
#' @param message Character string with progress message
#' @param step Current step number (optional)
#' @param total Total number of steps (optional)
#' @export
#' @examples
#' logProgress("Loading data files")
#' logProgress("Processing gene 1", step = 1, total = 1000)
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
#' and improve readability during execution.
#'
#' @param title Character string with section title
#' @export
#' @examples
#' logSection("Loading and Quality Control")
#' logSection("Pseudotime Analysis")
logSection <- function(title) {
  message("=== ", title, " ===")
}

# =============================================================================
# Interactive Selection Menus
# =============================================================================

#' Display interactive cell type selection menu
#'
#' Presents a numbered menu of available cell types for user selection.
#' Only functions in interactive R sessions. Provides clear error messages
#' for invalid selections or non-interactive contexts.
#'
#' @param cellTypes Character vector of available cell type names
#' @return Selected cell type as character string
#' @export
#' @examples
#' # Interactive usage:
#' cellTypes <- c("keratinocytes", "fibroblasts", "endothelial_cells")
#' selected <- showCellTypeMenu(cellTypes)
showCellTypeMenu <- function(cellTypes) {
  
  # Input validation
  if (is.null(cellTypes) || length(cellTypes) == 0) {
    stop("No valid cell types provided for selection")
  }
  
  # Check for interactive session
  if (!interactive()) {
    stop("Cell type selection requires an interactive R session.\n",
         "Available cell types: ", paste(cellTypes, collapse = ", "))
  }
  
  # Clear console and display menu
  clearConsole()
  cat("\n")
  
  # Present menu options
  selectionIndex <- menu(
    choices = cellTypes,
    title = "Select a cell type for analysis:"
  )
  
  # Handle user cancellation
  if (selectionIndex == 0) {
    stop("No cell type selected. Analysis cancelled.")
  }
  
  # Get selected cell type
  selectedCellType <- cellTypes[selectionIndex]
  
  # Clear console and confirm selection
  clearConsole()
  message("Selected cell type: ", selectedCellType)
  message("")
  
  return(selectedCellType)
}

#' Display interactive trajectory selection menu
#'
#' Presents a numbered menu of available pseudotime trajectories for user
#' selection. Only functions in interactive R sessions.
#'
#' @param trajectories Character vector of available trajectory identifiers
#' @return Selected trajectory as character string
#' @export
#' @examples
#' # Interactive usage:
#' trajectories <- c("Y_1", "Y_2", "Y_3")
#' selected <- showTrajectoryMenu(trajectories)
showTrajectoryMenu <- function(trajectories) {
  
  # Input validation
  if (is.null(trajectories) || length(trajectories) == 0) {
    stop("No valid pseudotime trajectories provided for selection")
  }
  
  # Check for interactive session
  if (!interactive()) {
    stop("Trajectory selection requires an interactive R session.\n",
         "Available trajectories: ", paste(trajectories, collapse = ", "))
  }
  
  # Clear console and display menu
  clearConsole()
  cat("\n")
  
  # Present menu options
  selectionIndex <- menu(
    choices = trajectories,
    title = "Select a pseudotime trajectory for analysis:"
  )
  
  # Handle user cancellation
  if (selectionIndex == 0) {
    stop("No trajectory selected. Analysis cancelled.")
  }
  
  # Get selected trajectory
  selectedTrajectory <- trajectories[selectionIndex]
  
  # Clear console and confirm selection
  clearConsole()
  message("Selected trajectory: ", selectedTrajectory)
  message("")
  
  return(selectedTrajectory)
}

# =============================================================================
# Confirmation Dialogs
# =============================================================================

#' Display confirmation dialog for potentially destructive operations
#'
#' Asks user to confirm before proceeding with operations that might
#' overwrite existing files or take a long time to complete.
#'
#' @param message Confirmation message to display
#' @param defaultYes Whether default response should be "yes" (default: FALSE)
#' @return Logical indicating user confirmation (TRUE = proceed, FALSE = cancel)
#' @export
#' @examples
#' if (confirmAction("This will overwrite existing results. Continue?")) {
#'   # Proceed with operation
#' }
confirmAction <- function(message, defaultYes = FALSE) {
  
  # In non-interactive mode, use default
  if (!interactive()) {
    message("Non-interactive mode: using default response (", 
            ifelse(defaultYes, "YES", "NO"), ")")
    return(defaultYes)
  }
  
  # Present confirmation dialog
  cat("\n")
  response <- menu(
    choices = c("Yes", "No"),
    title = paste0(message, "\nProceed?")
  )
  
  # Handle response
  if (response == 0) {
    # User cancelled menu
    return(FALSE)
  }
  
  return(response == 1)  # 1 = Yes, 2 = No
}

#' Display warning with user acknowledgment
#'
#' Shows a warning message and waits for user acknowledgment before
#' continuing script execution.
#'
#' @param message Warning message to display
#' @param waitForUser Whether to wait for user input (default: TRUE)
#' @export
#' @examples
#' showWarningDialog("Analysis may take several hours to complete")
showWarningDialog <- function(message, waitForUser = TRUE) {
  
  cat("\n")
  cat("WARNING: ", message, "\n")
  
  if (interactive() && waitForUser) {
    cat("Press [Enter] to continue or [Ctrl+C] to cancel...")
    readLines(n = 1)
  }
  
  cat("\n")
}

# =============================================================================
# Parameter Input Functions
# =============================================================================

#' Get numeric parameter from user input
#'
#' Prompts user for a numeric parameter with validation and default values.
#' Includes range checking and retry logic for invalid inputs.
#'
#' @param prompt Message to display for parameter request
#' @param default Default value if user provides no input
#' @param minValue Minimum allowed value (optional)
#' @param maxValue Maximum allowed value (optional)
#' @param maxRetries Maximum number of retry attempts (default: 3)
#' @return Numeric value from user input or default
#' @export
#' @examples
#' cores <- getNumericParameter("Number of cores to use", default = 1, minValue = 1, maxValue = 8)
getNumericParameter <- function(prompt, default = NULL, minValue = NULL, maxValue = NULL, maxRetries = 3) {
  
  # In non-interactive mode, use default
  if (!interactive()) {
    if (is.null(default)) {
      stop("Non-interactive mode requires default value for parameter: ", prompt)
    }
    message("Non-interactive mode: using default value (", default, ") for: ", prompt)
    return(default)
  }
  
  attempts <- 0
  
  while (attempts < maxRetries) {
    # Create prompt message
    promptMsg <- prompt
    if (!is.null(default)) {
      promptMsg <- paste0(promptMsg, " [default: ", default, "]")
    }
    if (!is.null(minValue) || !is.null(maxValue)) {
      rangeMsg <- ""
      if (!is.null(minValue) && !is.null(maxValue)) {
        rangeMsg <- paste0(" (", minValue, "-", maxValue, ")")
      } else if (!is.null(minValue)) {
        rangeMsg <- paste0(" (>= ", minValue, ")")
      } else {
        rangeMsg <- paste0(" (<= ", maxValue, ")")
      }
      promptMsg <- paste0(promptMsg, rangeMsg)
    }
    promptMsg <- paste0(promptMsg, ": ")
    
    # Get user input
    cat(promptMsg)
    userInput <- trimws(readLines(n = 1))
    
    # Use default if empty input
    if (userInput == "" && !is.null(default)) {
      return(default)
    }
    
    # Try to convert to numeric
    numericValue <- suppressWarnings(as.numeric(userInput))
    
    if (is.na(numericValue)) {
      cat("Invalid input. Please enter a numeric value.\n")
      attempts <- attempts + 1
      next
    }
    
    # Check range constraints
    if (!is.null(minValue) && numericValue < minValue) {
      cat("Value too small. Minimum allowed:", minValue, "\n")
      attempts <- attempts + 1
      next
    }
    
    if (!is.null(maxValue) && numericValue > maxValue) {
      cat("Value too large. Maximum allowed:", maxValue, "\n")
      attempts <- attempts + 1
      next
    }
    
    # Valid input
    return(numericValue)
  }
  
  # Max retries exceeded
  stop("Maximum retry attempts exceeded for parameter: ", prompt)
}

#' Get text parameter from user input
#'
#' Prompts user for a text parameter with validation and default values.
#' Includes basic validation and retry logic.
#'
#' @param prompt Message to display for parameter request
#' @param default Default value if user provides no input
#' @param allowEmpty Whether empty strings are acceptable (default: FALSE)
#' @param validOptions Optional vector of valid options to choose from
#' @param maxRetries Maximum number of retry attempts (default: 3)
#' @return Character string from user input or default
#' @export
#' @examples
#' method <- getTextParameter("Analysis method", default = "standard", 
#'                           validOptions = c("standard", "fast", "comprehensive"))
getTextParameter <- function(prompt, default = NULL, allowEmpty = FALSE, validOptions = NULL, maxRetries = 3) {
  
  # In non-interactive mode, use default
  if (!interactive()) {
    if (is.null(default)) {
      stop("Non-interactive mode requires default value for parameter: ", prompt)
    }
    message("Non-interactive mode: using default value (", default, ") for: ", prompt)
    return(default)
  }
  
  attempts <- 0
  
  while (attempts < maxRetries) {
    # Create prompt message
    promptMsg <- prompt
    if (!is.null(default)) {
      promptMsg <- paste0(promptMsg, " [default: ", default, "]")
    }
    if (!is.null(validOptions)) {
      promptMsg <- paste0(promptMsg, " (", paste(validOptions, collapse = ", "), ")")
    }
    promptMsg <- paste0(promptMsg, ": ")
    
    # Get user input
    cat(promptMsg)
    userInput <- trimws(readLines(n = 1))
    
    # Use default if empty input
    if (userInput == "" && !is.null(default)) {
      return(default)
    }
    
    # Check if empty input is allowed
    if (userInput == "" && !allowEmpty) {
      cat("Empty input not allowed. Please provide a value.\n")
      attempts <- attempts + 1
      next
    }
    
    # Check valid options
    if (!is.null(validOptions) && !userInput %in% validOptions) {
      cat("Invalid option. Valid choices are:", paste(validOptions, collapse = ", "), "\n")
      attempts <- attempts + 1
      next
    }
    
    # Valid input
    return(userInput)
  }
  
  # Max retries exceeded
  stop("Maximum retry attempts exceeded for parameter: ", prompt)
}

# =============================================================================
# Progress Bar Utilities
# =============================================================================

#' Create and manage simple text progress bar
#'
#' Creates a simple text-based progress bar for long-running operations.
#' Updates can be called incrementally to show progress.
#'
#' @param total Total number of operations
#' @param width Width of progress bar in characters (default: 50)
#' @return Environment containing progress bar functions
#' @export
#' @examples
#' pb <- createProgressBar(100)
#' for (i in 1:100) {
#'   # Do work here
#'   pb$update(i)
#'   Sys.sleep(0.1)
#' }
#' pb$finish()
createProgressBar <- function(total, width = 50) {
  
  current <- 0
  startTime <- Sys.time()
  
  # Progress bar environment
  pb <- new.env()
  
  # Update function
  pb$update <- function(value, message = NULL) {
    current <<- value
    
    # Calculate progress
    progress <- current / total
    filled <- floor(progress * width)
    empty <- width - filled
    
    # Create bar string
    bar <- paste0("[", 
                  paste(rep("=", filled), collapse = ""),
                  paste(rep(" ", empty), collapse = ""),
                  "]")
    
    # Calculate timing
    elapsed <- difftime(Sys.time(), startTime, units = "secs")
    if (current > 0) {
      eta <- (elapsed / current) * (total - current)
      etaStr <- sprintf("ETA: %s", format(.POSIXct(eta, tz = "UTC"), "%H:%M:%S"))  
    } else {
      etaStr <- "ETA: --:--:--"
    }
    
    # Format output
    output <- sprintf("\r%s %d/%d (%.1f%%) %s", 
                      bar, current, total, progress * 100, etaStr)
    
    if (!is.null(message)) {
      output <- paste0(output, " - ", message)
    }
    
    cat(output)
    flush.console()
  }
  
  # Finish function
  pb$finish <- function(message = "Complete") {
    pb$update(total, message)
    cat("\n")
  }
  
  return(pb)
}