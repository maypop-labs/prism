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
# Interactive Selection Menus
# =============================================================================

#' Display interactive cell type selection menu
#'
#' Presents a numbered menu of available cell types for user selection.
#' If only one cell type is available, automatically selects it.
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
  
  # Auto-select if only one option
  if (length(cellTypes) == 1) {
    clearConsole()
    message("Only one cell type available: ", cellTypes[1])
    message("Auto-selecting: ", cellTypes[1])
    message("")
    return(cellTypes[1])
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
#' selection. If only one trajectory is available, automatically selects it.
#' Only functions in interactive R sessions.
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
  
  # Auto-select if only one option
  if (length(trajectories) == 1) {
    clearConsole()
    message("Only one trajectory available: ", trajectories[1])
    message("Auto-selecting: ", trajectories[1])
    message("")
    return(trajectories[1])
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
