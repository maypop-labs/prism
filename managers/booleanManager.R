# -----------------------------------------------------------------------------
# Cell Type Menu
# -----------------------------------------------------------------------------

showCellTypeMenu <- function(cellTypes) {
  # Prompt the user only in an interactive R session
  if (is.null(cellTypes) || length(cellTypes) == 0) {
    stop("There are no valid cell types.")
  }
  if (interactive()) {
    cat("\014")
    cat("\n")
    selectionIndex <- menu(cellTypes, title = "Select a cell type:")
    
    if (selectionIndex == 0) {
      stop("No selection made. Exiting.")
    }
    cellType <- cellTypes[selectionIndex]
    cat("\014")
    message("You chose: ", cellType)
    message(" ")
  } else {
    stop("This script must be run interactively to choose a cell type.")
  }
  cellType
}

# -----------------------------------------------------------------------------
# Pseudotime Trajectory Menu
# -----------------------------------------------------------------------------

showTrajectoryMenu <- function(trajectories) {
  # Prompt the user only in an interactive R session
  if (is.null(trajectories) || length(trajectories) == 0) {
    stop("There are no valid pseudotime trajectories for this cell type.")
  }
  
  if (interactive()) {
    cat("\014")
    cat("\n")
    selectionIndex <- menu(trajectories, title = "Select a pseudotime trajectory:")
    
    if (selectionIndex == 0) {
      stop("No selection made. Exiting.")
    }
    traj <- trajectories[selectionIndex]
    cat("\014")
    message("You chose: ", traj)
    message(" ")
  } else {
    stop("This script must be run interactively to choose a cell type.")
  }
  traj
}
