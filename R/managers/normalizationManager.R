# =============================================================================
# normalizationManager.R
# Purpose: Donor normalization within cell types to prevent over-representation
# Dependencies: Seurat
# =============================================================================

# =============================================================================
# Core Normalization Function
# =============================================================================

#' Normalize donor representation within each cell type
#'
#' Downsamples cells within each (donor × cell type) stratum so that all donors
#' contribute equally to each cell type. Uses the minimum donor count as the
#' target, rejecting cell types where the minimum falls below a threshold.
#'
#' The function operates on cell-type-annotated Seurat objects and equalizes
#' donor representation separately within each cell type. This prevents
#' pseudotime inference and downstream analyses from being biased by donors
#' that contribute disproportionately many cells to specific cell types.
#'
#' @param seuratObj Seurat object with donor and cell type metadata
#' @param donorColumn Name of metadata column containing donor IDs (default: "donorID")
#' @param cellTypeColumn Name of metadata column containing cell types (default: "cellType")
#' @param floorCells Minimum acceptable cells per donor per cell type (default: 200)
#' @param numReplicates Number of independent random samplings (default: 1)
#' @param seedValue Random seed for reproducibility (default: 42)
#'
#' @return If numReplicates = 1: Single normalized Seurat object
#'         If numReplicates > 1: List of normalized Seurat objects
#'
#' @details
#' Target Calculation:
#' For each cell type, the target is the minimum donor count. Cell types where
#' this minimum is below floorCells are excluded from the normalized object.
#' This ensures perfect donor balance: all retained cell types have exactly
#' the same number of cells from each donor.
#'
#' Sampling Strategy:
#' - Find minimum donor count for each cell type
#' - If minimum >= floorCells: downsample all donors to match minimum
#' - If minimum < floorCells: exclude that cell type entirely
#' - Result: Perfect donor balance or complete exclusion
#'
#' Metadata Storage:
#' Normalization parameters are stored in the @misc slot of returned Seurat
#' objects for downstream validation.
#'
#' @export
normalizeCellTypeByDonor <- function(
  seuratObj,
  donorColumn = "donorID",
  cellTypeColumn = "cellType",
  floorCells = 200,
  numReplicates = 1,
  seedValue = 42
) {
  
  # --- Input Validation ---
  if (!donorColumn %in% colnames(seuratObj@meta.data)) {
    stop("Donor column '", donorColumn, "' not found in Seurat metadata")
  }
  
  if (!cellTypeColumn %in% colnames(seuratObj@meta.data)) {
    stop("Cell type column '", cellTypeColumn, "' not found in Seurat metadata")
  }
  
  if (floorCells < 1) {
    stop("floorCells must be >= 1")
  }
  
  if (numReplicates < 1) {
    stop("numReplicates must be >= 1")
  }
  
  # --- Calculate Target Cells Per Cell Type ---
  set.seed(seedValue)
  
  # Build contingency table: cell types × donors
  counts <- table(seuratObj@meta.data[[cellTypeColumn]], 
                  seuratObj@meta.data[[donorColumn]])
  
  # Calculate target as minimum donor count for each cell type
  nTarget <- apply(counts, 1, function(x) {
    nonzero <- x[x > 0]
    if (length(nonzero) == 0) return(0)
    
    # Use actual minimum donor count
    return(min(nonzero))
  })
  
  # Filter out cell types below floor threshold
  validCellTypes <- names(nTarget)[nTarget >= floorCells]
  excludedCellTypes <- names(nTarget)[nTarget < floorCells & nTarget > 0]
  
  if (length(excludedCellTypes) > 0) {
    message("\n=== Cell Types Excluded (Below Floor Threshold) ===")
    for (ct in excludedCellTypes) {
      donorCounts <- counts[ct, counts[ct, ] > 0]
      message("  ", ct, ": minimum ", nTarget[[ct]], " cells/donor (< ", 
              floorCells, " threshold)")
    }
  }
  
  if (length(validCellTypes) == 0) {
    stop("No cell types have sufficient cells for normalization. ",
         "Minimum required: ", floorCells, " cells per donor. ",
         "Consider lowering floorCells in config.yaml.")
  }
  
  # Keep only valid cell types
  nTarget <- nTarget[validCellTypes]
  
  # --- Log Normalization Strategy ---
  message("\n=== Donor Normalization Within Cell Types ===")
  message("Strategy: Use minimum donor count as target for perfect balance")
  message("Minimum threshold: ", floorCells, " cells per donor per cell type")
  message("Replicates: ", numReplicates)
  message("\nTargets per cell type:")
  
  for (ct in names(nTarget)) {
    donorCounts <- counts[ct, counts[ct, ] > 0]
    message("  ", ct, ": ", nTarget[[ct]], " cells/donor (from ",
            length(donorCounts), " donors with ", 
            paste(range(donorCounts), collapse = "-"), " cells)")
  }
  
  # --- Generate Replicates ---
  replicateList <- list()
  
  for (r in seq_len(numReplicates)) {
    if (numReplicates > 1) {
      set.seed(seedValue + r)
      message("\nGenerating replicate ", r, "/", numReplicates)
    }
    
    # Sample cells for each valid cell type only
    selectedCells <- unlist(lapply(names(nTarget), function(ct) {
      
      # Get unique donors in this cell type
      donors <- unique(seuratObj@meta.data[[donorColumn]][
        seuratObj@meta.data[[cellTypeColumn]] == ct
      ])
      
      target <- nTarget[[ct]]
      
      # Sample from each donor
      unlist(lapply(donors, function(dn) {
        theseCells <- rownames(subset(
          seuratObj@meta.data,
          seuratObj@meta.data[[cellTypeColumn]] == ct & 
          seuratObj@meta.data[[donorColumn]] == dn
        ))
        
        if (length(theseCells) > target) {
          # Downsample to target
          sample(theseCells, target)
        } else {
          # Keep all cells (no upsampling)
          theseCells
        }
      }))
    }))
    
    # Subset Seurat object to selected cells
    replicateObj <- subset(seuratObj, cells = selectedCells)
    
    # Store normalization metadata
    replicateObj@misc$normalization_meta <- list(
      nTarget = nTarget,
      seed = seedValue + r,
      floorCells = floorCells,
      totalCells = length(selectedCells),
      originalCells = ncol(seuratObj),
      retentionRate = length(selectedCells) / ncol(seuratObj),
      excludedCellTypes = if(length(excludedCellTypes) > 0) excludedCellTypes else NULL,
      donorColumn = donorColumn,
      cellTypeColumn = cellTypeColumn
    )
    
    replicateList[[r]] <- replicateObj
  }
  
  # --- Report Results ---
  message("\n=== Normalization Complete ===")
  message("Original cells: ", ncol(seuratObj))
  message("Retained cells: ", ncol(replicateList[[1]]))
  message("Retention rate: ", 
          round(100 * ncol(replicateList[[1]]) / ncol(seuratObj), 1), "%")
  
  # Validate donor balance
  if (numReplicates == 1) {
    validateDonorBalance(replicateList[[1]], donorColumn, cellTypeColumn, nTarget)
  }
  
  # --- Return Results ---
  if (numReplicates == 1) {
    return(replicateList[[1]])
  } else {
    names(replicateList) <- paste0("rep", seq_len(numReplicates))
    return(replicateList)
  }
}

# =============================================================================
# Validation and Reporting
# =============================================================================

#' Validate donor balance after normalization
#'
#' Checks that each donor's representation per cell type is exactly equal
#' to the target (perfect balance). Reports any discrepancies.
#'
#' @param seuratObj Normalized Seurat object
#' @param donorColumn Name of donor metadata column
#' @param cellTypeColumn Name of cell type metadata column
#' @param nTarget Named vector of target counts per cell type
#'
#' @return Invisible TRUE if validation passes
#' @keywords internal
validateDonorBalance <- function(seuratObj, donorColumn, cellTypeColumn, nTarget) {
  
  message("\n=== Validating Donor Balance ===")
  
  # Rebuild contingency table after normalization
  countsAfter <- table(seuratObj@meta.data[[cellTypeColumn]], 
                       seuratObj@meta.data[[donorColumn]])
  
  allPerfect <- TRUE
  
  for (ct in names(nTarget)) {
    donorCounts <- countsAfter[ct, countsAfter[ct, ] > 0]
    target <- nTarget[[ct]]
    
    # Check if all donors have exactly the target count
    if (length(unique(donorCounts)) == 1 && unique(donorCounts) == target) {
      message("  \u2713 ", ct, ": Perfect balance (", target, " cells per donor)")
    } else {
      allPerfect <- FALSE
      message("  \u2717 ", ct, ": Imperfect balance")
      message("    Target: ", target, ", Actual: ", 
              paste(donorCounts, collapse = ", "))
    }
  }
  
  if (allPerfect) {
    message("\n\u2713 All cell types have perfect donor balance")
  } else {
    warning("Some cell types do not have perfect balance. ",
            "This should not happen with the minimum-based approach.")
  }
  
  invisible(TRUE)
}

#' Extract normalization metadata from Seurat object
#'
#' Retrieves normalization parameters stored in @misc slot for reporting
#' or validation purposes.
#'
#' @param seuratObj Normalized Seurat object
#'
#' @return List of normalization metadata, or NULL if not normalized
#' @export
getNormalizationMetadata <- function(seuratObj) {
  if (!"normalization_meta" %in% names(seuratObj@misc)) {
    warning("No normalization metadata found in Seurat object")
    return(NULL)
  }
  
  return(seuratObj@misc$normalization_meta)
}

#' Report donor statistics per cell type
#'
#' Generates a summary table of donor representation across cell types,
#' useful for understanding dataset composition before and after normalization.
#'
#' @param seuratObj Seurat object with donor and cell type metadata
#' @param donorColumn Name of donor metadata column (default: "donorID")
#' @param cellTypeColumn Name of cell type metadata column (default: "cellType")
#' @param showDetailed Logical indicating whether to print detailed per-donor counts (default: TRUE)
#'
#' @return Data frame with donor statistics per cell type
#' @export
reportDonorStats <- function(
  seuratObj,
  donorColumn = "donorID",
  cellTypeColumn = "cellType",
  showDetailed = TRUE
) {
  
  # Build contingency table
  counts <- table(seuratObj@meta.data[[cellTypeColumn]], 
                  seuratObj@meta.data[[donorColumn]])
  
  # Calculate statistics per cell type
  stats <- data.frame(
    cellType = rownames(counts),
    totalCells = rowSums(counts),
    numDonors = apply(counts, 1, function(x) sum(x > 0)),
    minCells = apply(counts, 1, function(x) min(x[x > 0])),
    maxCells = apply(counts, 1, function(x) max(x[x > 0])),
    medianCells = apply(counts, 1, function(x) median(x[x > 0])),
    stringsAsFactors = FALSE
  )
  
  # Sort by total cells descending
  stats <- stats[order(stats$totalCells, decreasing = TRUE), ]
  rownames(stats) <- NULL
  
  # Print detailed per-donor counts if requested
  if (showDetailed) {
    message("\n=== Per-Donor Cell Counts by Cell Type ===")
    for (ct in stats$cellType) {
      donorCounts <- counts[ct, ]
      donorCounts <- donorCounts[donorCounts > 0]  # Only show donors with cells
      message("\n", ct, " (total: ", stats$totalCells[stats$cellType == ct], ")")
      for (dn in names(donorCounts)) {
        message("  ", dn, ": ", donorCounts[dn], " cells")
      }
    }
    message("\n")
  }
  
  return(stats)
}
