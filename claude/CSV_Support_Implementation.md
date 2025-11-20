# CSV Support Implementation for PRISM

## Date: November 2025

## Overview

PRISM now supports CSV format count matrices in addition to the existing HDF5 (.h5) and CellRanger directory formats. This enables analysis of datasets from sources like the Human Cell Atlas that provide pre-processed marker expression files in CSV format.

## Implementation Details

### Files Modified

1. **`pathManager.R`** - Added CSV detection and reading capabilities
2. **`01 - Loading, Aggregation, and Quality Control.R`** - Added CSV loading logic

### Key Changes

#### 1. `buildCellrangerPaths()` Function

**Location**: `R/managers/pathManager.R`

**Changes**:
- Added CSV file pattern matching using `Sys.glob()`
- Implemented priority order: `.h5` > directory > `.csv`
- Updated error messages to show all three expected path formats

**Pattern Matching**: Uses `[donorID]_marker_expr*.csv` or `[donorID]_output.csv` pattern to match files

#### 2. New `readCsvMatrix()` Function

**Location**: `R/managers/pathManager.R`

**Purpose**: Reads CSV matrix files and converts to sparse matrix format

**Details**:
- Expects genes as rows, cell barcodes as columns
- First row contains cell barcodes
- First column contains gene names
- Automatically converts to sparse matrix format for memory efficiency

**Function Signature**:
```r
readCsvMatrix <- function(csvPath)
```

**Returns**: Sparse matrix compatible with `CreateSeuratObject()`

#### 3. Script 01 Loading Logic

**Location**: `R/01 - Loading, Aggregation, and Quality Control.R`

**Changes**:
- Restructured file type detection to use nested if-else with extension checking
- Added `.csv` detection via `grepl("\\.csv$", paths$cellranger[i])`
- Calls `readCsvMatrix()` for CSV files

**Loading Flow**:
1. Check if path is a directory → use `Read10X()`
2. Check if path is a file:
   - If `.h5` extension → use `Read10X_h5()`
   - If `.csv` extension → use `readCsvMatrix()`
   - Otherwise → error for unsupported file type

## Usage

### File Naming Conventions

The CSV files should follow one of these naming patterns:

1. **Standard pattern**: `[donorID]_output.csv`
   - Example: `26yo_output.csv`, `56yo_output.csv`, `67yo_output.csv`

2. **Extended pattern**: `[donorID]_marker_expr*.csv`
   - Example: `26yo_marker_exprLOH1.csv`, `56yo_marker_exprHA2.csv`

### CSV File Format Requirements

- **Structure**: Genes in rows, cells in columns
- **First row**: Cell barcodes (e.g., `AAACCCAAG...`)
- **First column**: Gene names/symbols
- **Content**: Integer count values (sparse matrix, mostly zeros)

### config.yaml Setup

No changes needed to existing config.yaml structure. The system will automatically detect and use CSV files:

```yaml
rootPath: "E:/datasets/omics/testis/"

donorIDs:
  - "26yo"
  - "56yo"
  - "67yo"
```

## Benefits

1. **Flexibility**: Supports mixed formats within same project
2. **Backward Compatible**: Existing `.h5` and directory workflows unchanged
3. **Memory Efficient**: Converts CSV to sparse matrix immediately
4. **Automated**: No manual file conversion required
5. **Clear Error Messages**: Shows all expected path formats when files not found

## Testing

### Testis Dataset

**Location**: `E:/datasets/omics/testis/cellranger_counts/`

**Files**:
- `26yo_output.csv`
- `56yo_output.csv`
- `67yo_output.csv`

**Format**: Human Cell Atlas marker expression matrices (genes × cells)

### Expected Workflow

1. User places CSV files in `cellranger_counts/` folder
2. Files named to match `donorIDs` in `config.yaml`
3. Run Script 01 normally
4. System automatically detects CSV format and uses `readCsvMatrix()`
5. Downstream analysis proceeds identically to `.h5`/directory workflows

## Technical Notes

### Priority Order Rationale

1. **`.h5` first**: Most efficient format, fastest to read, smallest file size
2. **Directory second**: Traditional CellRanger output, well-tested
3. **`.csv` last**: Fallback for special datasets, less efficient but flexible

### Sparse Matrix Conversion

The `readCsvMatrix()` function immediately converts the dense CSV matrix to sparse format because:
- scRNA-seq data is typically ~95% zeros
- Sparse format dramatically reduces memory usage
- Seurat expects sparse matrices internally
- No information loss in conversion

### Pattern Matching

Uses `Sys.glob()` for flexible pattern matching:
- Handles wildcards in filenames
- Returns all matching files
- Takes first match if multiple files found
- Useful for datasets with complex suffixes (e.g., LOH1, HA2)

## Future Enhancements

Possible future improvements:
1. Support for gzipped CSV files (`.csv.gz`)
2. Progress bar for large CSV reading
3. Validation of CSV structure before loading
4. Support for alternative delimiters (TSV, etc.)
5. Batch conversion of CSV to HDF5 for repeated analyses

## Conclusion

PRISM now seamlessly handles three input formats, making it compatible with a wider range of publicly available datasets while maintaining the same user experience and downstream analysis pipeline.
