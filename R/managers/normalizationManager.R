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
#' Downsamples cells within each (donor × cell type) stratum so that no donor
#' dominates any given cell type. Uses a soft quantile-based target rather than
#' strict minimum to preserve statistical power while reducing donor bias.
#'
#' The function operates on cell-type-annotated Seurat objects and equalizes
#' donor representation separately within each cell type. This prevents
#' pseudotime inference and downstream analyses from being biased by donors
#' that contribute disproportionately many cells to specific cell types.
#'
#' @param seuratObj Seurat object with donor and cell type metadata
#' @param donorColumn Name of metadata column containing donor IDs (default: "donorID")
#' @param cellTypeColumn Name of metadata column containing cell types (default: "cellType")
#' @param floorCells Minimum target cells per donor per cell type (default: 200)
#' @param capCells Maximum target cells per donor per cell type (default: 1000)
#' @param quantileCut Quantile for soft target calculation (default: 0.2)
#' @param numReplicates Number of independent random samplings (default: 1)
#' @param seedValue Random seed for reproducibility (default: 42)
#'
#' @return If numReplicates = 1: Single normalized Seurat object
#'         If numReplicates > 1: List of normalized Seurat objects
#'
#' @details
#' Target Calculation:
#' For each cell type, the target is the 20th percentile of non-zero donor
#' counts, bounded by floorCells and capCells. This avoids both:
#' - Over-downsampling (using strict minimum)
#' - Extreme donors dominating (using mean or median)
#'
#' Sampling Strategy:
#' - Donors with more cells than target: randomly sample target cells
#' - Donors with fewer cells than target: keep all cells (no upsampling)
#'
#' Metadata Storage:
#' Normalization parameters are stored in the @misc slot of returned Seurat
#' objects for downstream validation and leave-one-donor-out (LODO) analysis.
#'
#' @export
normalizeCellTypeByDonor <- function(
  seuratObj,
  donorColumn = "donorID",
  cellTypeColumn = "cellType",
  floorCells = 200,
  capCells = 1000,
  quantileCut = 0.2,
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
  
  if (capCells < floorCells) {
    stop("capCells must be >= floorCells")
  }
  
  if (quantileCut < 0 || quantileCut > 1) {
    stop("quantileCut must be between 0 and 1")
  }
  
  if (numReplicates < 1) {
    stop("numReplicates must be >= 1")
  }
  
  # --- Calculate Target Cells Per Cell Type ---
  set.seed(seedValue)
  
  # Build contingency table: cell types × donors
  counts <- table(seuratObj@meta.data[[cellTypeColumn]], 
                  seuratObj@meta.data[[donorColumn]])
  
  # Calculate soft target for each cell type
  nTarget <- apply(counts, 1, function(x) {
    nonzero <- x[x > 0]
    if (length(nonzero) == 0) return(0)
    
    # Use quantileCut percentile, bounded by floor and cap
    target <- quantile(nonzero, quantileCut)
    target <- max(min(target, capCells), floorCells)
    return(round(target))
  })
  
  # --- Log Normalization Strategy ---
  message("\n=== Donor Normalization Within Cell Types ===")
  message("Strategy: ", quantileCut * 100, "th percentile of donor counts")
  message("Bounds: [", floorCells, ", ", capCells, "] cells per donor per cell type")
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
    
    # Sample cells for each cell type
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
      capCells = capCells,
      quantileCut = quantileCut,
      totalCells = length(selectedCells),
      originalCells = ncol(seuratObj),
      retentionRate = length(selectedCells) / ncol(seuratObj),
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
#' Checks that each donor's representation per cell type is within expected
#' range after normalization. Reports any cell types where balance is poor.
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
  
  imbalanced <- character()
  
  for (ct in names(nTarget)) {
    donorCounts <- countsAfter[ct, countsAfter[ct, ] > 0]
    target <- nTarget[[ct]]
    
    # Check if any donor deviates >20% from target
    deviations <- abs(donorCounts - target) / target
    
    if (any(deviations > 0.2)) {
      imbalanced <- c(imbalanced, ct)
      message("  WARNING: ", ct, " has donors deviating >20% from target")
      message("    Target: ", target, ", Range: ", 
              paste(range(donorCounts), collapse = "-"))
    }
  }
  
  if (length(imbalanced) == 0) {
    message("  All cell types within ±20% of target")
  } else {
    message("\nCell types with imbalance: ", paste(imbalanced, collapse = ", "))
    message("This may indicate insufficient cells for balanced sampling")
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
#'
#' @return Data frame with donor statistics per cell type
#' @export
reportDonorStats <- function(
  seuratObj,
  donorColumn = "donorID",
  cellTypeColumn = "cellType"
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
  
  return(stats)
}
