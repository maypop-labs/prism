# PRISM: Pseudotime-Resolved Identification of Switch-gene Mechanisms

![PRISM Logo](Logo.png)

PRISM is a comprehensive bioinformatics pipeline for analyzing cellular aging using single-cell RNA sequencing data. The pipeline integrates pseudotime analysis, gene regulatory network construction, Boolean attractor modeling, and perturbation analysis to discover genes that could shift cellular states toward more youthful configurations.

## Overview

PRISM identifies potential therapeutic targets for cellular rejuvenation by:
1. Inferring aging trajectories from single-cell RNA-seq data
2. Detecting genes with switch-like behavior during aging
3. Constructing gene regulatory networks using SCENIC/GENIE3
4. Building Boolean network models with attractor dynamics
5. Testing perturbations to identify rejuvenation targets

## Quick Start

### 1. Environment Setup

**Run this script once before any other analysis:**

```r
source("R/setup_environment.R")
```

This script installs all necessary R package dependencies, including Bioconductor packages, CRAN packages, and development versions from GitHub.

### 2. Data Preparation

Single-cell RNA sequencing data should be:
- **Format**: Gene-by-cell transcript counts
- **Location**: `{rootPath}/counts/` folder
- **Accepted formats**:
  - CellRanger output folders (containing `matrix.mtx`, `features.tsv`, `barcodes.tsv`)
  - `.h5` files (HDF5 format)
  - `.csv` files (gene-by-cell matrices)

**File Naming Convention:**
Files or folders must be named using the `donorID` values from `config.yaml`, followed by `_output`:
- Example: `26yo_output.h5`
- Example: `56yo_output` (folder with CellRanger files)
- Example: `67yo_output.csv`

### 3. Configuration

Edit `R/config.yaml` to set up your experiment. **The three essential parameters are:**

```yaml
rootPath: "E:/datasets/omics/your_tissue/"  # Base directory for your project
ages: [26, 56, 67]                           # Ages of your donors (in years)
donorIDs: ["26yo", "56yo", "67yo"]          # Identifiers matching your data files
```

Additional important parameters:
- `rcisTargetPath`: Directory containing RcisTarget databases for SCENIC
- `cores`: Number of CPU cores to use for parallel processing
- `saveResults`: Whether to save intermediate and final results (recommended: `true`)

See [YAML Documentation Comment Standard.md](YAML%20Documentation%20Comment%20Standard.md) for detailed documentation of all configuration parameters.

### 4. Run the Pipeline

Execute the numbered R scripts in order from the `R/` directory:

```r
# From R console or RStudio
setwd("E:/bin/prism/R")

source("01 - Loading, Aggregation, and Quality Control.R")
source("02 - Cell Type Identification.R")
source("03 - Pseudotime.R")
# ... and so on
```

Or run individual scripts as needed during development.

## Pipeline Scripts

The pipeline consists of 11 numbered scripts that should be run sequentially:

| Script | Purpose |
|--------|---------|
| **01 - Loading, Aggregation, and Quality Control** | Load scRNA-seq data, detect doublets, perform QC filtering, merge donors, and integrate with Seurat |
| **02 - Cell Type Identification** | Annotate cell types with SingleR, normalize donor representation, create cell-type-specific objects |
| **03 - Pseudotime** | Infer trajectories with Monocle3, validate age-pseudotime correlations, identify aging-related branches |
| **04 - GeneSwitches** | Identify genes with switch-like expression transitions along pseudotime using logistic regression |
| **05 - GRN** | Construct gene regulatory networks using SCENIC (GENIE3 + RcisTarget motif analysis) |
| **06 - Boolean Regulation** | Transform GRN into Boolean logic rules, optimize rules based on expression data |
| **07 - BoolNet** | Create executable Boolean network, identify attractor states and basin statistics |
| **08 - Attractor Analysis** | Compute aging scores for attractors, test single and double gene perturbations for rejuvenation targets |
| **09 - Graphs** | Generate publication-ready visualizations of results (attractor landscapes, perturbation rankings) |
| **10 - Switch Gene Regulatory Report** | Analyze regulatory relationships for switch genes, create detailed reports |
| **11 - Switch Gene Perturbation Analysis** | Test switch gene perturbations specifically, complementing full-network analysis |

See [PRISM - Project Overview.md](PRISM%20-%20Project%20Overview.md) in the `claude/` folder for detailed descriptions of each script.

## Output Structure

The pipeline creates a `results/` folder at `{rootPath}` with the following structure:

```
results/
└── celltype_name/
   ├── rds/           # R data objects (.rds files)
   ├── tsv/           # Tab-separated value tables
   ├── txt/           # Text reports and logs
   ├── plots/         # Visualization plots (.png files)
   ├── rds/           # R data files
   ├── tsv/           # Tab-separated value files (tables)
   ├── txt/           # Output logs
   ├── plots/         # plots
   └── scenic/        # Intermediate files (SCENIC checkpoints)
```

## Key Features

### Donor Normalization
The pipeline implements a critical donor normalization strategy that achieves **perfect balance** of donor representation within each cell type. This prevents donor-specific biases in downstream analyses while preserving biological variation. See `normalizeDonors` and `normFloorCells` parameters in `config.yaml`.

### Trajectory Validation
Pseudotime trajectories are rigorously validated by requiring minimum correlations (r ≥ 0.3) between pseudotime and chronological age. This ensures trajectories represent genuine aging processes rather than noise.

### Checkpoint Recovery
Long-running analyses (especially SCENIC/GRN construction) include checkpoint systems to enable resumption from interruptions. Intermediate results are stored in `int/` directories.

### Comprehensive Documentation
All configuration parameters are thoroughly documented using the [YAML Documentation Comment Standard](YAML%20Documentation%20Comment%20Standard.md), providing practical guidance with suggested values and rationales.

## Requirements

### Software
- R (≥ 4.0.0)
- RStudio (recommended)
- Required R packages installed via `setup_environment.R`

### Data
- Single-cell RNA-seq data in CellRanger, HDF5, or CSV format
- RcisTarget databases for SCENIC analysis (download from https://resources.aertslab.org/cistarget/)
  - `hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather`
  - `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather`

### Computational Resources
- **Minimum**: 8GB RAM, 4 CPU cores
- **Recommended**: 32GB RAM, 8 CPU cores
- **Large datasets**: 64-128GB RAM, 16 CPU cores
- **Disk space**: 50-200GB depending on dataset size

## Typical Runtime

Runtimes vary greatly with dataset size and system resources:
- **Scripts 01-02**: 30 minutes - 2 hours
- **Script 03**: 20 minutes - 1 hour per cell type
- **Script 04**: 10-30 minutes per trajectory
- **Script 05 (SCENIC)**: 14-48 hours per trajectory (longest step)
- **Scripts 06-08**: 2-8 hours per trajectory
- **Scripts 09-11**: 10-30 minutes per trajectory

Total pipeline runtime: **1-3 days per cell type** depending on complexity.

## Troubleshooting

### Common Issues

**"Cannot find cellranger_counts folder"**
- Ensure `rootPath` in config.yaml points to the directory containing `counts/` subfolder
- Check that files are named correctly: `{donorID}_output.h5` or `{donorID}_output/`

**"No trajectories passed correlation threshold"**
- Consider lowering `pseudotimeMinAgeCorrelation` (minimum 0.3)
- Check if age range in your data is sufficient (>20 years recommended)
- Verify cell type has enough cells across all age groups

**"SCENIC analysis failed or incomplete"**
- Ensure `rcisTargetPath` points to directory with .feather database files
- Check that you have sufficient RAM (SCENIC is memory-intensive)
- Use checkpoint recovery: set `grnFromScratch: false` to resume

**"Out of memory errors"**
- Reduce `seuratNumberOfFeatures` or `singleRNumberOfFeatures`
- Lower `cores` to reduce parallel memory usage
- Increase `maxMemoryGB` parameter (but stay below system RAM)

### Getting Help

For questions, issues, or contributions:
1. Check the detailed script documentation in `claude/` folder
2. Review configuration parameter documentation in config.yaml
3. Consult the [R Style Guide](R%20Style%20Guide.md) for code conventions

## Citation

If you use PRISM in your research, please cite:
- **PRISM pipeline**: [Publication pending]
- **Supporting methods**:
  - Seurat: Hao et al. (2021) Cell
  - Monocle3: Cao et al. (2019) Nature
  - GeneSwitches: Psarras et al. (2021) Bioinformatics
  - SCENIC: Aibar et al. (2017) Nature Methods
  - BoolNet: Müssel et al. (2010) Bioinformatics

**Last Updated**: November 2025  
**Version**: 1.0
