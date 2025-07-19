# PRISM
## Pseudotime Reversion via In Silico Modeling

This pipeline performs pseudotime-based regulatory analysis of single-cell RNA-seq data, inferring dynamic gene regulatory networks and assessing candidate aging-reversal targets through Boolean modeling and attractor perturbation.

## Overview
The workflow consists of 15 sequential R scripts, each performing a defined step:

### Step-by-Step Breakdown

1. **01 - Loading, Aggregation, and Quality Control.R**
   - Loads raw 10X datasets from multiple individuals
   - Merges and filters Seurat objects based on mitochondrial content and feature count
   - Normalizes, clusters, and computes PCA + UMAP embeddings

2. **02 - Cell Type Identification.R**
   - Uses `SingleR` and `HumanPrimaryCellAtlasData()` to annotate cell types
   - Subsets top 3 abundant types for downstream analysis

3. **03 - Pseudotime Core.R**
   - Converts Seurat object to Monocle3 `CellDataSet`
   - Learns graph and pseudotime ordering from a specified root cell type
   - Extracts trajectory branches that correlate with age (r > 0.6)

4. **04 - DEGs.R**
   - Performs pseudotime-based graph tests to identify differentially expressed genes (DEGs)

5. **05 - Smoothing.R**
   - Smooths expression matrix across pseudotime using a moving average

6. **06 - GeneSwitches.R**
   - Applies logistic regression to find switch-like expression transitions
   - Filters top 1000 switch-like DEGs

7. **07 - GRN - Part 01.R**
   - Constructs a co-expression and motif-based GRN using SCENIC

8. **08 - GRN - Part 02.R**
   - Filters and signs SCENIC GRN edges by Spearman correlation
   - Removes terminal nodes and merges strongly connected components

9. **09 - Boolean Regulation.R**
   - Infers Boolean logic rules for each gene using binarized data and regulatory structure

10. **10 - BoolNet.R**
    - Converts Boolean rules to a BoolNet model
    - Computes all attractors using stochastic sampling

11. **11 - Attractor Analysis - Part 01.R**
    - Projects attractors onto pseudotime-derived young–old axis
    - Scores each attractor for aging state

12. **12 - Attractor Analysis - Part 02.R**
    - Computes entropy and stability for each attractor under perturbation
    - Combines age, stability, and basin size into a composite score

13. **13 - Attractor Analysis - Part 03.R**
    - Systematically perturbs genes (single and pairwise)
    - Measures change in aging score to identify reversion targets

14. **14 - Conclusion.R**
    - Compiles perturbation results into a clean summary TSV

15. **15 - Graphs.R**
    - Generates publication-quality figures for trajectory, GRN, attractors, and perturbation scores

## Output Files
- `*.rds`: Intermediate R objects (e.g., Seurat, Monocle3, GRN graphs)
- `*.tsv`: Summary statistics and perturbation scores
- `*.png`: Visualization figures

## Helper Scripts
- `functions.R`: Defines all shared utility functions (entropy, Boolean rule parsing, attractor decoding, etc.)

## Flowchart
(see generated diagram)

## Requirements
- R 4.2+
- Packages: Seurat, monocle3, BoolNet, SCENIC, GeneSwitches, igraph, tidyverse

## Contact
For usage questions or contributions, contact: [Your Name / Lab / Email]

