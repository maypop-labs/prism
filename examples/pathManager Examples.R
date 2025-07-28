# =============================================================================
# Example Usage Functions
# =============================================================================

#' Example: How a script with no dependencies would use this system
exampleScript01Usage <- function() {
  # For script 01 - Loading, Aggregation, and Quality Control
  pathInfo <- initializeInteractivePaths(needsCellType = FALSE, needsTrajectory = FALSE)
  config <- pathInfo$config
  paths <- pathInfo$paths
  
  ensureProjectDirectories(paths)
  
  # Now use paths like:
  # mergedSeurat <- readRDS(paths$static$mergedSeurat)
  # ggsave(paths$static$pcaAllCells, plot)
  
  return(list(config = config, paths = paths))
}

#' Example: How a script requiring cell type would use this system
exampleScript03Usage <- function() {
  # For script 03 - Pseudotime Core
  pathInfo <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = FALSE)
  config <- pathInfo$config
  paths <- pathInfo$paths
  cellType <- pathInfo$cellType
  
  ensureProjectDirectories(paths)
  
  # Now use paths like:
  # seuratObj <- readRDS(paths$cellType$seuratObject)
  # saveRDS(retainedLeaves, paths$cellType$retainedTrajectories)
  
  return(list(config = config, paths = paths, cellType = cellType))
}

#' Example: How a script requiring both cell type and trajectory would use this system
exampleScript09Usage <- function() {
  # For script 09 - Boolean Regulation
  pathInfo <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
  config <- pathInfo$config
  paths <- pathInfo$paths
  cellType <- pathInfo$cellType
  trajectory <- pathInfo$trajectory
  
  ensureProjectDirectories(paths)
  
  # Now use paths like:
  # cds <- load_monocle_objects(directory_path = paths$trajectory$monocle3SmoothedGeneSwitches)
  # boolRules <- readRDS(paths$trajectory$booleanRules)
  
  return(list(config = config, paths = paths, cellType = cellType, trajectory = trajectory))
}