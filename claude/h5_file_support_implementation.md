# H5 File Support Implementation

**Date**: November 17, 2025  
**Purpose**: Add support for loading 10X Genomics data from `.h5` files in addition to existing folder structure support

## Summary

Modified PRISM pipeline to support both CellRanger output formats:
1. **Traditional format**: Folder structure with matrix files (`filtered_feature_bc_matrix/`)
2. **HDF5 format**: Single `.h5` files (e.g., `23yo_output.h5`)

## Changes Made

### 1. Modified `pathManager.R` - `buildCellrangerPaths()` function

**Location**: `E:/bin/prism/R/managers/pathManager.R`

**Change**: Rewrote function to auto-detect file format for each donor

```r
buildCellrangerPaths <- function(config) {
  basePath <- paste0(config$rootPath, "cellranger_counts/")
  paths <- character(length(config$donorIDs))
  
  for (i in seq_along(config$donorIDs)) {
    # Check for .h5 file first (newer format, more efficient)
    h5Path <- paste0(basePath, config$donorIDs[i], "_output.h5")
    dirPath <- paste0(basePath, config$donorIDs[i], "_output/outs/filtered_feature_bc_matrix/")
    
    if (file.exists(h5Path)) {
      paths[i] <- h5Path
    } else if (dir.exists(dirPath)) {
      paths[i] <- dirPath
    } else {
      stop("No valid CellRanger data found for ", config$donorIDs[i], 
           ". Expected either:\n  ", h5Path, "\n  or: ", dirPath)
    }
  }
  
  return(paths)
}
```

**Logic**:
- Prioritizes `.h5` files if both formats exist (newer, more efficient)
- Falls back to directory structure if `.h5` not found
- Provides clear error message showing expected paths if neither found

### 2. Modified Script 01 - Data Loading Logic

**Location**: `E:/bin/prism/R/01 - Loading, Aggregation, and Quality Control.R`

**Change**: Added file type detection and appropriate loading function selection

```r
# Detect file type and use appropriate loading function
if (dir.exists(paths$cellranger[i])) {
  # Load from CellRanger directory structure
  counts <- Read10X(data.dir = paths$cellranger[i])
} else if (file.exists(paths$cellranger[i]) && grepl("\\.h5$", paths$cellranger[i])) {
  # Load from HDF5 file
  counts <- Read10X_h5(filename = paths$cellranger[i])
} else {
  stop("Invalid CellRanger path: ", paths$cellranger[i])
}
```

**Logic**:
- Uses `Read10X()` for directory-based data
- Uses `Read10X_h5()` for `.h5` files
- Validates path exists before attempting to load

## Benefits

1. **Backward Compatibility**: Existing analyses using folder structure continue to work
2. **Future Compatibility**: Supports newer HDF5 format from 10X Genomics
3. **Flexibility**: Can mix formats within same project (some donors as folders, others as .h5)
4. **Efficiency**: HDF5 files are more compact and faster to load
5. **Clear Errors**: Helpful error messages if data not found

## Testing Recommendations

1. Verify `.h5` files exist in `E:/datasets/omics/ovary/cellranger_counts/`
2. Check file naming matches pattern: `{donorID}_output.h5` (e.g., `23yo_output.h5`)
3. Run Script 01 to verify loading works correctly
4. Check console output for confirmation of which format is being loaded

## Expected File Locations

For the ovary dataset with 8 donors (23yo through 54yo):

```
E:/datasets/omics/ovary/cellranger_counts/
├── 23yo_output.h5
├── 27yo_output.h5
├── 28yo_output.h5
├── 29yo_output.h5
├── 49yo_output.h5
├── 51yo_output.h5
├── 52yo_output.h5
└── 54yo_output.h5
```

## Notes

- The `.h5` format is HDF5-based and requires the `hdf5r` package (dependency of Seurat)
- Function `Read10X_h5()` is built into Seurat and handles the HDF5 format
- No changes needed to config.yaml - the detection is automatic
- No changes needed to downstream scripts - they work with Seurat objects regardless of input format
