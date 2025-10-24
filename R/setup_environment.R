# =============================================================================
# setup_environment
#
# Loads and installs a list of necessary R packages.
#
# =============================================================================

options(warn = -1)

# --- Functions ---

script_path <- function() {
  ## 1.  Was the script run with Rscript ... --file=path ?
  cmd_args   <- commandArgs(trailingOnly = FALSE)
  file_arg   <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1)
    return(normalizePath(sub("^--file=", "", file_arg)))
  
  ## 2.  Was it sourced with source("path/to/file.R") ?
  if (!is.null(sys.frame(1)$ofile))
    return(normalizePath(sys.frame(1)$ofile))
  
  ## 3.  Are we inside RStudio with the file open?
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- rstudioapi::getSourceEditorContext()
    if (nzchar(ctx$path))
      return(normalizePath(ctx$path))
  }
  
  ## 4.  Give up politely
  stop("Cannot figure out where this script lives.")
}

# --- Script Execution ---

script_dir    <- paste0(dirname(script_path()), "/")
user_dir      <- paste0(gsub("\\\\", "/", normalizePath("~")), "/")
lib_path      <- paste0(user_dir, "R/", "win-library/", R.version$major, ".", R.version$minor, "/")
renviron_path <- paste0(user_dir, ".Renviron")
renviron_line <- paste0('R_LIBS_USER="', lib_path, '"')
rprofile_path <- paste0(user_dir, ".Rprofile")

options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Create a .Renviron file ---
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
writeLines(renviron_line, renviron_path)
.libPaths(c(lib_path, .libPaths()))

# --- Create an .Rprofile file ---

rprofile_contents <- character()
if (file.exists(rprofile_path)) {
  rprofile_contents <- readLines(rprofile_path, warn = FALSE)
}
libpath_line <- sprintf('.libPaths(c("%s", .libPaths()))', lib_path)

# Only add it if not already present
if (!any(grepl("libPaths", rprofile_contents))) {
  writeLines(c(rprofile_contents, "", "# Set user library path", libpath_line), rprofile_path)
  message("Added .libPaths() override to your .Rprofile at: ", rprofile_path)
} else {
  message(".Rprofile already contains a .libPaths() setting. Skipping edit.")
}

# --- Read in the list of required packages ---

pkgs <- readLines(paste0(script_dir, "required_packages.txt"))
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
} else { message("All packages are already installed.") }


# --- Manually Install BP Cells ---

install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))

# --- Read in the list of bioconductor packages ---

pkgs <- readLines(paste0(script_dir, "bioconductor_packages.txt"))
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) {
  BiocManager::install(new_pkgs)
} else { message("All devtool packages are already installed.") }

# --- Read in the list of devtool packages ---

pkgs <- readLines(paste0(script_dir, "devtools_packages.txt"))
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) {
  devtools::install_github(new_pkgs)
} else { message("All devtool packages are already installed.") }


message("Done!")
