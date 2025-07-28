buildPath <- function(type, filename = NULL, cellType = NULL, trajectory = NULL) {
  basePath <- paste0(config$rootPath, "results/", type, "/")
  if (!is.null(filename)) {
    if (!is.null(cellType) && !is.null(trajectory)) {
      paste0(basePath, cellType, "_", trajectory, "_", filename)
    } else if (!is.null(cellType)) {
      paste0(basePath, cellType, "_", filename)
    } else {
      paste0(basePath, filename)
    }
  } else {
    basePath
  }
}

initializeScript <- function() {
  options(warn = -1)
  config <- yaml::read_yaml("config.yaml")
  options(Seurat.object.assay.version = config$SeuratAssay)
  registerDoParallel(cores = config$cores)
  loadRequiredLibraries()
  return(config)
}


loadRequiredLibraries <- function() {
  standardLibs <- c(
    "BoolNet",
    "celldex",
    "data.table",
    "doParallel",
    "dplyr",
    "foreach",
    "GeneSwitches",
    "GENIE3",
    "ggplot2",
    "ggraph",
    "gmp",
    "igraph",
    "Matrix",
    "mgcv",
    "monocle3",
    "progress",
    "progressr",
    "purrr",
    "readr",
    "SCENIC",
    "Seurat",
    "SeuratData",
    "SeuratDisk",
    "SingleCellExperiment",
    "SingleR",
    "tidyr",
    "tidyverse",
    "tools",
    "yaml"
  )
  lapply(standardLibs, library, character.only = TRUE)
  invisible(TRUE)
}

timeIt <- function(expr, name = "task") {
  start <- Sys.time()
  result <- force(expr)
  end   <- Sys.time()
  dur   <- difftime(end, start, units = "secs")
  cat(sprintf("%s finished in %d min %.1f s\n",
              name, as.integer(dur) %/% 60, as.numeric(dur) %% 60))
  invisible(result)
}