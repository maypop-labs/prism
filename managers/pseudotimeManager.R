# =============================================================================
# pseudotimeManager.R
# Purpose: Functions for analyzing Monocle3 pseudotime trajectories
# =============================================================================

calculateTrajectoryCorrelation <- function(subCds) {
  # Extract data
  pt <- pseudotime(subCds)
  ages <- colData(subCds)$age
  donors <- colData(subCds)$donor_id
  
  # Clean data
  valid_idx <- !is.na(pt) & is.finite(pt) & !is.na(ages)
  
  if (sum(valid_idx) < 10) {
    return(list(correlation = NA, p_value = NA, method = "insufficient_data"))
  }
  
  pt_clean <- pt[valid_idx]
  ages_clean <- ages[valid_idx]
  donors_clean <- donors[valid_idx]
  n_donors <- length(unique(donors_clean))
  
  # Primary test: Spearman correlation (tests monotonic relationship)
  spearman_test <- cor.test(pt_clean, ages_clean, method = "spearman")
  
  # Secondary test: Mixed-effects model (if enough donors)
  lme_result <- NULL
  if (n_donors >= 4) {
    tryCatch({
      df <- data.frame(
        pseudotime = pt_clean, 
        age = ages_clean, 
        donor_id = as.factor(donors_clean)
      )
      lme_model <- lmer(pseudotime ~ age + (1|donor_id), data = df)
      lme_summary <- summary(lme_model)
      
      lme_result <- list(
        age_coefficient = lme_summary$coefficients["age", "Estimate"],
        age_pvalue = lme_summary$coefficients["age", "Pr(>|t|)"],
        converged = TRUE
      )
    }, error = function(e) { 
      lme_result <- list(age_coefficient = NA, age_pvalue = NA, converged = FALSE)
    })
  }
  
  # Choose primary result based on available data
  if (!is.null(lme_result) && lme_result$converged) {
    primary_correlation = lme_result$age_coefficient  # Effect size from mixed model
    primary_p_value = lme_result$age_pvalue
    method_used = "mixed_effects"
  } else {
    primary_correlation = spearman_test$estimate
    primary_p_value = spearman_test$p.value  
    method_used = "spearman_correlation"
  }
  
  return(list(
    correlation = primary_correlation,
    p_value = primary_p_value,
    method = method_used,
    spearman_rho = spearman_test$estimate,
    spearman_p = spearman_test$p.value,
    lme_result = lme_result,
    n_cells = sum(valid_idx),
    n_donors = n_donors
  ))
}

convertSeuratToCDS <- function(seuratObj) {
  rawCounts <- LayerData(seuratObj, layer = "counts")
  genes <- rownames(rawCounts)
  cellMeta <- seuratObj@meta.data
  geneMeta <- data.frame(gene_id = genes, gene_short_name = genes, row.names = genes)
  new_cell_data_set(rawCounts, cellMeta, geneMeta)
}

getRootCells <- function(cds, rootAge) {
  meta <- as.data.frame(colData(cds))
  rownames(meta)[meta$age == rootAge]
}

runPseudotime <- function(cds, rootCells) {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "log", scaling = TRUE, verbose = config$verbose)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose = config$verbose)
  cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden", verbose = config$verbose)
  cds <- learn_graph(cds)
  order_cells(cds, reduction_method = "UMAP", root_cells = rootCells)
}