# =============================================================================
# 05 - GRN Debug
#
# Debug version of 05 - GRN.R to identify the matrix dimension issue
# =============================================================================

# --- Initialization ---
source("managers/attractorManager.R")
source("managers/booleanManager.R")
source("managers/grnManager.R")
source("managers/pathManager.R")
source("managers/pseudotimeManager.R")
source("managers/setupManager.R")
source("managers/uiManager.R")

config     <- initializeScript()
pathInfo   <- initializeInteractivePaths(needsCellType = TRUE, needsTrajectory = TRUE)
paths      <- pathInfo$paths
cellType   <- pathInfo$cellType
trajectory <- pathInfo$trajectory
ctPaths    <- getCellTypeFilePaths(paths$base, cellType)
ptPaths    <- getTrajectoryFilePaths(paths$base, cellType, trajectory)
ensureProjectDirectories(paths)
clearConsole()

# --- Load pseudotime trajectory ---
if (!dir.exists(ptPaths$monocle3GeneSwitches)) stop("Monocle3 object directory not found: ", ptPaths$monocle3GeneSwitches)
if (config$verbose) { message("Loading Monocle3 object from: ", ptPaths$monocle3GeneSwitches) }
cds <- load_monocle_objects(directory_path = ptPaths$monocle3GeneSwitches)

# DEBUG: Check CDS object
cat("DEBUG: CDS dimensions:", dim(cds), "\n")
cat("DEBUG: Available assays:", assayNames(cds), "\n")

# --- Load Switch Genes ---
if (!file.exists(ptPaths$geneSwitches)) stop("Switch genes RDS file not found: ", ptPaths$geneSwitches)
if (config$verbose) { message("Loading switch genes from: ", ptPaths$geneSwitches) }
switchGenes <- readRDS(ptPaths$geneSwitches)

# DEBUG: Check switch genes
cat("DEBUG: Switch genes class:", class(switchGenes), "\n")
cat("DEBUG: Switch genes length/nrow:", ifelse(is.data.frame(switchGenes), nrow(switchGenes), length(switchGenes)), "\n")
if (is.data.frame(switchGenes)) {
  cat("DEBUG: Switch genes column names:", colnames(switchGenes), "\n")
  cat("DEBUG: First few rownames:", head(rownames(switchGenes)), "\n")
} else {
  cat("DEBUG: Switch genes structure:", str(switchGenes), "\n")
}

# --- Filter and Normalize Expression Matrix ---
switchGeneNames <- rownames(switchGenes)
cat("DEBUG: Number of switch gene names:", length(switchGeneNames), "\n")
cat("DEBUG: First few switch gene names:", head(switchGeneNames), "\n")

exprMat <- assay(cds, "counts")
cat("DEBUG: Initial exprMat dimensions:", dim(exprMat), "\n")
cat("DEBUG: exprMat class:", class(exprMat), "\n")

keepGenes <- rowSums(exprMat > 1) >= 10
cat("DEBUG: Number of genes passing filter (>1 counts in >=10 cells):", sum(keepGenes), "\n")

exprMat <- exprMat[keepGenes, ]
cat("DEBUG: exprMat dimensions after filtering:", dim(exprMat), "\n")

cellInfo <- data.frame(row.names = colnames(exprMat))

# --- Initialize SCENIC ---
data(defaultDbNames)
species  <- config$grnScenicSpecies
baseName <- paste0("motifAnnotations_", species)
# Ordered list of possible dataset names shipped in different SCENIC releases
candidates <- c(baseName, paste0(baseName, "_v10"), paste0(baseName, "_v9"), "motifAnnotations")
loaded <- FALSE
for (ds in candidates) {
  if (ds %in% data(package = "RcisTarget")$results[, "Item"]) {
    data(list = ds, package = "RcisTarget", envir = environment())
    assign(baseName, get(ds), envir = globalenv())
    loaded <- TRUE
    break
  }
}
if (!loaded) { stop("No motif annotation object found for '", species, "'. Check your RcisTarget installation or download the annotation manually.")  }

if (config$verbose) { message("Initializing SCENIC") }
scenicOptions <- initializeScenic(
  org    = config$grnScenicSpecies,
  dbDir  = config$rcisTargetPath,
  dbs    = config$grnScenicDBs,
  nCores = config$cores
)

# --- Filter to switch genes and TFs ---
allTFs     <- getDbTfs(scenicOptions)
cat("DEBUG: Number of TFs from database:", length(allTFs), "\n")

unionGenes <- union(switchGeneNames, allTFs)
cat("DEBUG: Number of union genes (switch + TFs):", length(unionGenes), "\n")

intersectGenes <- intersect(rownames(exprMat), unionGenes)
cat("DEBUG: Number of genes in intersection:", length(intersectGenes), "\n")

exprMat <- exprMat[intersectGenes, ]
cat("DEBUG: Final exprMat dimensions before SCENIC:", dim(exprMat), "\n")
cat("DEBUG: Final exprMat class:", class(exprMat), "\n")

# Check if exprMat is valid for rowSums
if (length(dim(exprMat)) < 2) {
  stop("ERROR: exprMat has insufficient dimensions: ", paste(dim(exprMat), collapse="x"))
}

if (nrow(exprMat) == 0) {
  stop("ERROR: No genes remain after filtering")
}

if (ncol(exprMat) == 0) {
  stop("ERROR: No cells remain after filtering")
}

# --- Run SCENIC Pipeline ---
if (config$verbose) { message("Filtering genes and calculating co-expression") }

# DEBUG: Check right before the problematic function call
cat("DEBUG: About to call geneFiltering with exprMat dimensions:", dim(exprMat), "\n")
cat("DEBUG: exprMat is matrix:", is.matrix(exprMat), "\n")
cat("DEBUG: exprMat is array:", is.array(exprMat), "\n")

# Try the problematic function call
tryCatch({
  filteredGenes <- geneFiltering(exprMat, scenicOptions)
  cat("DEBUG: Successfully filtered genes. Number remaining:", length(filteredGenes), "\n")
}, error = function(e) {
  cat("ERROR in geneFiltering:", e$message, "\n")
  cat("DEBUG: exprMat summary:\n")
  print(summary(exprMat))
  stop("Stopping due to geneFiltering error")
})

message("Debug completed successfully!")