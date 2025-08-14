# =============================================================================
# 12 - Conclusion
#
# Consolidate single- and double-gene perturbation results and prepare for output.
# Writes combined results to TSV files for downstream interpretation.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(readr)

# --- Parameters ---
cellType       <- "Keratinocytes"
cellTrajectory <- "Y_447"

# --- Input Files ---
s0File <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_single_targets_0.rds")
s1File <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_single_targets_1.rds")
d0File <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_double_targets_KD_KD.rds")
d1File <- paste0("E:/datasets/omics/skin/results/rds/", cellType, "_", cellTrajectory, "_final_double_targets_OE_OE.rds")

# --- Output Files ---
singleOutFile <- paste0("E:/datasets/omics/skin/results/tsv/", cellType, "_single_gene_reversion_summary.tsv")
doubleOutFile <- paste0("E:/datasets/omics/skin/results/tsv/", cellType, "_double_gene_reversion_summary.tsv")

# --- Load Data ---
s0 <- if (file.exists(s0File)) readRDS(s0File) else NULL
s1 <- if (file.exists(s1File)) readRDS(s1File) else NULL
d0 <- if (file.exists(d0File)) readRDS(d0File) else NULL
d1 <- if (file.exists(d1File)) readRDS(d1File) else NULL

# --- Combine Single-Gene Perturbations ---
singleCombined <- bind_rows(
  mutate(s0, Mode = "Knockdown"),
  mutate(s1, Mode = "Overexpression")
) %>% pivot_wider(
  names_from = Mode,
  values_from = AgingScore,
  names_prefix = "AgingScore_"
) %>% arrange(pmin(AgingScore_Knockdown, AgingScore_Overexpression, na.rm = TRUE))

# --- Combine Double-Gene Perturbations ---
doubleCombined <- bind_rows(
  mutate(d0, Mode = "Knockdown"),
  mutate(d1, Mode = "Overexpression")
) %>% pivot_wider(
  names_from = Mode,
  values_from = AgingScore,
  names_prefix = "AgingScore_"
) %>% arrange(pmin(AgingScore_Knockdown, AgingScore_Overexpression, na.rm = TRUE))

# --- Save Outputs ---
if (!is.null(singleCombined)) write_tsv(singleCombined, singleOutFile)
if (!is.null(doubleCombined)) write_tsv(doubleCombined, doubleOutFile)

# --- Optional Display ---
if (!is.null(singleCombined)) print(head(singleCombined, 10))
if (!is.null(doubleCombined)) print(head(doubleCombined, 10))

message("Done!")
