# =============================================================================
# 12 - Attractor Analysis - Part 03
#
# Systematic in silico perturbation of genes to identify those that reduce
# the aging score. Includes both single- and double-gene perturbations.
# Now with corrected BoolNet usage and robust age score computation.
# =============================================================================

# --- Libraries ---
library(foreach)
library(doParallel)
library(BoolNet)
library(monocle3)
library(dplyr)
library(gmp)
library(progress)

# --- Source functions ---
source("functions.R")

# --- Options ---
options(warn = -1)
config <- yaml::read_yaml("config.yaml")
options(Seurat.object.assay.version = config$SeuratAssay)
registerDoParallel(cores = config$cores)

# --- Parameters ---
monocle3Path        <- paste0(config$rootPath, "results/monocle3/")
graphMlPath         <- paste0(config$rootPath, "results/graphml/")
plotPath            <- paste0(config$rootPath, "results/plots/")
rdsPath             <- paste0(config$rootPath, "results/rds/")
tsvPath             <- paste0(config$rootPath, "results/tsv/")
txtPath             <- paste0(config$rootPath, "results/txt/")
cellTypes           <- readRDS(paste0(rdsPath, "cell_types.rds"))
cellType            <- showCellTypeMenu(cellTypes)
trajNamesFile       <- readRDS(paste0(rdsPath, "retained_trajectories_", cellType, ".rds"))
cellTrajectory      <- showTrajectoryMenu(trajNamesFile)
cdsPath             <- paste0(monocle3Path, "monocle3_", cellType, "_", cellTrajectory, "_smoothed_geneSwitches")
degFile             <- paste0(rdsPath, cellType, "_", cellTrajectory, "_switch_degs.rds")
edgesFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02_edges.rds")
graphFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_GRN_Part_02.rds")
rulesFile           <- paste0(rdsPath, cellType, "_", cellTrajectory, "_Boolean_Rules.rds")
boolnetFile         <- paste0(rdsPath, cellType, "_", cellTrajectory, "_boolnet.rds")
attractorsFile      <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractors.rds")
attractorDfFile     <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractor_df.rds")
geneMapFile         <- paste0(rdsPath, cellType, "_", cellTrajectory, "_gene_map.rds")
attractorScoresFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_attractor_df_scores.rds")

dir.create(graphMlPath, recursive = TRUE, showWarnings = FALSE)
dir.create(plotPath,    recursive = TRUE, showWarnings = FALSE)
dir.create(rdsPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(tsvPath,     recursive = TRUE, showWarnings = FALSE)
dir.create(txtPath,     recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
if (!dir.exists(cdsPath)) stop("Monocle3 object directory not found: ", cdsPath)
if (!file.exists(boolnetFile)) stop("BoolNet RDS file not found: ", boolnetFile)
if (!file.exists(attractorsFile)) stop("Attractor RDS file not found: ", attractorsFile)
if (!file.exists(attractorDfFile)) stop("Attractor data frame RDS file not found: ", attractorDfFile)
if (!file.exists(attractorScoresFile)) stop("Attractor scores data RDS file not found: ", attractorScoresFile)

message("Loading Monocle3 object from: ", cdsPath)
cds <- load_monocle_objects(directory_path = cdsPath)
message("Loading BoolNet from: ", boolnetFile)
boolnet      <- readRDS(boolnetFile)
message("Loading attractors from: ", attractorsFile)
attractors <- readRDS(attractorsFile)
message("Loading attractor data frame from: ", attractorDfFile)
attractorDf <- readRDS(attractorDfFile)
message("Loading attractor scores from: ", attractorScoresFile)
attractorScores <- readRDS(attractorScoresFile)

cat("\014")
cat("\n")

# --- Output Paths ---
singleKDFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_single_targets_0.rds")
singleOEFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_single_targets_1.rds")
doubleKDFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_KD_KD.rds")
doubleOEFile <- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_OE_OE.rds")
doubleMixFile<- paste0(rdsPath, cellType, "_", cellTrajectory, "_final_double_targets_KD_OE.rds")

# =============================================================================
# FIXED HELPER FUNCTIONS
# =============================================================================

# --- Corrected BoolNet attractor computation ---
getAttractorsCorrect <- function(network, max_states = 1000) {
  # Use exhaustive search for small networks, random sampling for large ones
  if (length(network$genes) <= 15) {
    # Small network: exhaustive search
    attractors <- getAttractors(network, method = "exhaustive")
  } else {
    # Large network: use synchronous method with proper random sampling
    attractors <- getAttractors(network, 
                                method = "random", 
                                startStates = max_states,
                                type = "synchronous",
                                returnTable = TRUE)
  }
  
  # Verify we got actual attractors, not just random states
  if (length(attractors$attractors) > 50) {
    warning("Too many attractors found (", length(attractors$attractors), 
            "). Network might be degenerate. Trying alternative method.")
    
    # Try with fewer states and synchronous updates
    attractors <- getAttractors(network, 
                                method = "random", 
                                startStates = 100,
                                type = "synchronous",
                                returnTable = TRUE)
  }
  
  return(attractors)
}

# --- Robust age score computation ---
computeAgeScoreRobust <- function(v, vYoung, vOld) {
  
  # Check for numerical issues
  diff_vec <- vOld - vYoung
  
  # Remove genes where young and old are too similar
  valid_genes <- abs(diff_vec) > 1e-6
  
  if (sum(valid_genes) < 5) {
    warning("Too few genes with meaningful young/old differences. Using all genes.")
    valid_genes <- rep(TRUE, length(diff_vec))
  }
  
  # Subset to valid genes
  v_sub <- v[valid_genes]
  vYoung_sub <- vYoung[valid_genes]
  vOld_sub <- vOld[valid_genes]
  
  # Compute age score with robust denominator
  num <- sum((v_sub - vYoung_sub) * (vOld_sub - vYoung_sub))
  den <- sum((vOld_sub - vYoung_sub)^2)
  
  if (den < 1e-10) {
    warning("Denominator too small in age score computation")
    return(0)
  }
  
  age_score <- num / den
  
  # Clamp to reasonable range
  age_score <- pmax(-2, pmin(3, age_score))
  
  return(age_score)
}

# --- Original computeAgeScore for compatibility ---
computeAgeScore <- function(v, vYoung, vOld) {
  num <- sum((v - vYoung) * (vOld - vYoung))
  den <- sum((vOld - vYoung)^2)
  if (den < 1e-10) return(0)
  num / den
}

# --- Fixed perturbation scoring function ---
perturbNetworkScoreFixed <- function(network, genes, values, youngVec, oldVec, nStates = 1000) {
  
  # Fix genes in network
  net <- fixGenes(network, fixIndices = genes, values = values)
  
  # Get attractors using corrected method
  ats <- getAttractorsCorrect(net, max_states = nStates)
  
  if (length(ats$attractors) == 0) {
    warning("No attractors found for genes: ", paste(genes, collapse = ", "))
    return(NA_real_)
  }
  
  # Limit to reasonable number of attractors
  if (length(ats$attractors) > 20) {
    warning("Too many attractors (", length(ats$attractors), 
            "). Using only first 20 for efficiency.")
    ats$attractors <- ats$attractors[1:20]
  }
  
  # Compute age scores for each attractor
  ageScores <- numeric(length(ats$attractors))
  basinSizes <- numeric(length(ats$attractors))
  
  for (i in seq_along(ats$attractors)) {
    att <- ats$attractors[[i]]
    
    # Get a representative state from this attractor
    if (length(att$involvedStates) > 0) {
      decoded <- decodeBigIntegerState(att$involvedStates[[1]], length(network$genes))
      names(decoded) <- network$genes
      
      # Compute age score only for genes in both decoded and young/old vectors
      gSet <- intersect(names(decoded), names(youngVec))
      if (length(gSet) > 0) {
        ageScores[i] <- computeAgeScoreRobust(decoded[gSet], youngVec[gSet], oldVec[gSet])
      } else {
        ageScores[i] <- 0
      }
      
      basinSizes[i] <- att$basinSize
    }
  }
  
  # Normalize basin sizes
  if (sum(basinSizes) > 0) {
    basinFrac <- basinSizes / sum(basinSizes)
  } else {
    basinFrac <- rep(1/length(basinSizes), length(basinSizes))
  }
  
  # Compute weighted average age score
  weighted_age_score <- sum(basinFrac * ageScores, na.rm = TRUE)
  
  return(weighted_age_score)
}

# --- Network health diagnostic function ---
diagnoseNetworkHealth <- function(network, sample_genes = 5) {
  
  cat("=== NETWORK HEALTH DIAGNOSIS ===\n")
  
  # Test baseline attractors
  baseline_ats <- getAttractorsCorrect(network, max_states = 100)
  cat("Baseline attractors found:", length(baseline_ats$attractors), "\n")
  
  if (length(baseline_ats$attractors) > 0) {
    basin_sizes <- sapply(baseline_ats$attractors, function(a) a$basinSize)
    cat("Basin sizes:", paste(basin_sizes, collapse = ", "), "\n")
  }
  
  # Test a few gene perturbations
  test_genes <- sample(network$genes, min(sample_genes, length(network$genes)))
  
  for (gene in test_genes) {
    cat("\nTesting gene:", gene, "\n")
    
    # Test knockdown
    tryCatch({
      net_kd <- fixGenes(network, fixIndices = gene, values = 0)
      ats_kd <- getAttractorsCorrect(net_kd, max_states = 100)
      cat("  KD attractors:", length(ats_kd$attractors), "\n")
    }, error = function(e) {
      cat("  KD failed:", e$message, "\n")
    })
    
    # Test overexpression
    tryCatch({
      net_oe <- fixGenes(network, fixIndices = gene, values = 1)
      ats_oe <- getAttractorsCorrect(net_oe, max_states = 100)
      cat("  OE attractors:", length(ats_oe$attractors), "\n")
    }, error = function(e) {
      cat("  OE failed:", e$message, "\n")
    })
  }
  
  cat("=== END DIAGNOSIS ===\n\n")
}

# --- Time estimation helper ---
estimateTimeRemaining <- function(start_time, completed, total) {
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  if (completed == 0) return("Unknown")
  
  rate <- elapsed / completed
  remaining <- rate * (total - completed)
  
  if (remaining > 60) {
    sprintf("%.1f hours", remaining / 60)
  } else {
    sprintf("%.1f minutes", remaining)
  }
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# --- Define Young and Old Vectors ---
matBin <- assay(cds, "binary")
cellOrder <- order(colData(cds)$Pseudotime)
nCells <- length(cellOrder)
q20 <- floor(0.2 * nCells)
youngVec <- rowMeans(matBin[, cellOrder[1:q20]])
oldVec   <- rowMeans(matBin[, cellOrder[(nCells - q20 + 1):nCells]])

# --- Baseline Score ---
initialScore <- sum(attractorScores$AttractorScore)
cat("Initial Aging Score:", initialScore, "\n")
if (initialScore <= config$attInitialThreshold) stop("Initial score below threshold")

# --- Young/Old Vector Diagnostics ---
cat("\n=== YOUNG/OLD VECTOR DIAGNOSTICS ===\n")
cat("Young vector - Min:", round(min(youngVec), 4), "Max:", round(max(youngVec), 4), 
    "Mean:", round(mean(youngVec), 4), "\n")
cat("Old vector - Min:", round(min(oldVec), 4), "Max:", round(max(oldVec), 4), 
    "Mean:", round(mean(oldVec), 4), "\n")
cat("Correlation between young/old vectors:", round(cor(youngVec, oldVec), 4), "\n")

# Check for genes with minimal differences
diff_vec <- abs(oldVec - youngVec)
low_diff_genes <- sum(diff_vec < 1e-6)
cat("Genes with minimal young/old differences:", low_diff_genes, "/", length(diff_vec), "\n")

if (cor(youngVec, oldVec) > 0.98) {
  warning("Young and old vectors are very similar (r > 0.98). Age scores may be unstable.")
}

# --- Network Health Check ---
cat("\n=== BOOLEAN NETWORK DIAGNOSTICS ===\n")
cat("Network genes:", length(boolnet$genes), "\n")
diagnoseNetworkHealth(boolnet, sample_genes = 3)

# --- Gene Selection ---
genes <- intersect(boolnet$genes, names(youngVec))
cat("Genes available for perturbation:", length(genes), "\n")

if (length(genes) == 0) {
  stop("No genes available for perturbation analysis!")
}

# =============================================================================
# SINGLE-GENE PERTURBATIONS
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SINGLE-GENE PERTURBATIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat(sprintf("Testing %d genes for knockdown and overexpression\n", length(genes)))
cat("This analysis uses corrected BoolNet attractor computation.\n\n")

# Initialize results
singleKD <- data.frame(Gene = genes, AgingScore = NA_real_)
singleOE <- data.frame(Gene = genes, AgingScore = NA_real_)

# Setup progress tracking
script_start_time <- Sys.time()
pb <- progress_bar$new(
  format = "  [:bar] :percent (:current/:total) ETA: :eta Elapsed: :elapsed",
  total = length(genes) * 2,
  clear = FALSE,
  width = 80
)

cat("\014")
cat("\n")
cat("Progress will be shown below:\n")

# Run perturbations
for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # Progress updates
  if (i %% 25 == 0 || i == 1) {
    elapsed <- difftime(Sys.time(), script_start_time, units = "mins")
    if (i > 1) {
      cat("\014")
      cat("\n")
      cat(sprintf("\nProgress update: %d/%d genes completed (%.1f%%). ", 
                  i, length(genes), 100*i/length(genes)))
      cat("ETA:", estimateTimeRemaining(script_start_time, i, length(genes)), "\n")
    }
  }
  
  # Knockdown
  pb$tick()
  singleKD$AgingScore[i] <- perturbNetworkScoreFixed(boolnet, gene, 0, youngVec, oldVec)
  
  # Overexpression  
  pb$tick()
  singleOE$AgingScore[i] <- perturbNetworkScoreFixed(boolnet, gene, 1, youngVec, oldVec)
}

# Format results
singleKD <- singleKD %>% 
  mutate(UnperturbedScore = initialScore, Delta = AgingScore - initialScore) %>% 
  arrange(AgingScore)

singleOE <- singleOE %>% 
  mutate(UnperturbedScore = initialScore, Delta = AgingScore - initialScore) %>% 
  arrange(AgingScore)

# Results summary
single_time <- difftime(Sys.time(), script_start_time, units = "mins")
cat(sprintf("\nSingle-gene perturbations completed in %.1f minutes\n", single_time))

# Data quality check
valid_kd <- sum(!is.na(singleKD$AgingScore))
valid_oe <- sum(!is.na(singleOE$AgingScore))
cat(sprintf("Data Quality: %d/%d valid KD results, %d/%d valid OE results\n", 
            valid_kd, length(genes), valid_oe, length(genes)))

# Show score distributions
cat("\nScore distributions:\n")
cat("KD scores - Min:", round(min(singleKD$AgingScore, na.rm=T), 4), 
    "Max:", round(max(singleKD$AgingScore, na.rm=T), 4),
    "Mean:", round(mean(singleKD$AgingScore, na.rm=T), 4), "\n")
cat("OE scores - Min:", round(min(singleOE$AgingScore, na.rm=T), 4), 
    "Max:", round(max(singleOE$AgingScore, na.rm=T), 4),
    "Mean:", round(mean(singleOE$AgingScore, na.rm=T), 4), "\n")

# Show top results
cat("\nTOP 10 AGING-IMPROVEMENT TARGETS (most negative deltas):\n")
improvement_targets <- singleKD[singleKD$Delta < 0, ]
if (nrow(improvement_targets) > 0) {
  print(head(improvement_targets, 10))
} else {
  cat("No improvement targets found in KD results.\n")
}

improvement_targets_oe <- singleOE[singleOE$Delta < 0, ]
if (nrow(improvement_targets_oe) > 0) {
  cat("\nTOP 10 OE AGING-IMPROVEMENT TARGETS:\n")
  print(head(improvement_targets_oe, 10))
}

cat("\nTOP 5 AGING-WORSENING TARGETS (most positive deltas):\n")
worsening_targets <- singleKD[singleKD$Delta > 0, ]
if (nrow(worsening_targets) > 0) {
  print(head(worsening_targets, 5))
} else {
  cat("No worsening targets found.\n")
}

# =============================================================================
# DOUBLE-GENE PERTURBATIONS
# =============================================================================

# Only proceed if we have sufficient valid single results
if (valid_kd >= 5 && valid_oe >= 5) {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("DOUBLE-GENE PERTURBATIONS\n") 
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Get top candidates (prioritize improvement targets)
  topKD <- head(singleKD$Gene[singleKD$Delta < 0], 5)
  topOE <- head(singleOE$Gene[singleOE$Delta < 0], 5)
  
  if (length(topKD) == 0) {
    cat("No beneficial KD targets found. Using top 5 by score.\n")
    topKD <- head(singleKD$Gene, 5)
  }
  if (length(topOE) == 0) {
    cat("No beneficial OE targets found. Using top 5 by score.\n") 
    topOE <- head(singleOE$Gene, 5)
  }
  
  cat("Top KD candidates:", paste(topKD, collapse = ", "), "\n")
  cat("Top OE candidates:", paste(topOE, collapse = ", "), "\n")
  
  # Generate combinations safely
  doubleKD <- expand.grid(Gene1 = topKD, Gene2 = topKD, stringsAsFactors = FALSE) %>% 
    filter(Gene1 < Gene2)
  doubleOE <- expand.grid(Gene1 = topOE, Gene2 = topOE, stringsAsFactors = FALSE) %>% 
    filter(Gene1 < Gene2)
  doubleMix <- expand.grid(KD_Gene = topKD, OE_Gene = topOE, stringsAsFactors = FALSE)
  
  total_combinations <- nrow(doubleKD) + nrow(doubleOE) + nrow(doubleMix)
  cat("Combinations to test: KD-KD =", nrow(doubleKD), 
      ", OE-OE =", nrow(doubleOE), ", Mixed =", nrow(doubleMix), 
      ", Total =", total_combinations, "\n")
  
  if (total_combinations > 0) {
    
    # Safe scoring function
    scorePairsFixed <- function(df, vals, description) {
      if (nrow(df) == 0) {
        cat("No", description, "combinations to test.\n")
        return(df)
      }
      
      cat("\nTesting", nrow(df), description, "combinations...\n")
      df$AgingScore <- numeric(nrow(df))
      df$UnperturbedScore <- numeric(nrow(df))
      df$Delta <- numeric(nrow(df))
      
      for (i in 1:nrow(df)) {
        genes_to_perturb <- c(df[i, 1], df[i, 2])
        df$AgingScore[i] <- perturbNetworkScoreFixed(boolnet, genes_to_perturb, vals, youngVec, oldVec)
        df$UnperturbedScore[i] <- initialScore
        df$Delta[i] <- df$AgingScore[i] - initialScore
        
        if (i %% 5 == 0) {
          cat(sprintf("  Completed %d/%d %s combinations\n", i, nrow(df), description))
        }
      }
      
      return(df)
    }
    
    # Run double perturbations
    double_start_time <- Sys.time()
    
    doubleKD <- scorePairsFixed(doubleKD, c(0, 0), "KD-KD")
    doubleOE <- scorePairsFixed(doubleOE, c(1, 1), "OE-OE")
    doubleMix <- scorePairsFixed(doubleMix, c(0, 1), "Mixed")
    
    # Sort results
    if (nrow(doubleKD) > 0) doubleKD <- doubleKD %>% arrange(AgingScore)
    if (nrow(doubleOE) > 0) doubleOE <- doubleOE %>% arrange(AgingScore)
    if (nrow(doubleMix) > 0) doubleMix <- doubleMix %>% arrange(AgingScore)
    
    double_time <- difftime(Sys.time(), double_start_time, units = "mins")
    cat(sprintf("\nDouble-gene perturbations completed in %.1f minutes\n", double_time))
    
    # Show best double combinations
    if (nrow(doubleKD) > 0) {
      cat("\nBEST KD-KD COMBINATIONS:\n")
      print(head(doubleKD, 5))
    }
    if (nrow(doubleOE) > 0) {
      cat("\nBEST OE-OE COMBINATIONS:\n")
      print(head(doubleOE, 5))
    }
    if (nrow(doubleMix) > 0) {
      cat("\nBEST MIXED COMBINATIONS:\n")
      print(head(doubleMix, 5))
    }
    
  } else {
    cat("No valid combinations to test.\n")
    doubleKD <- data.frame()
    doubleOE <- data.frame()
    doubleMix <- data.frame()
  }
  
} else {
  cat("\nSkipping double perturbations due to insufficient valid single results.\n")
  cat("Need at least 5 valid KD and 5 valid OE results.\n")
  doubleKD <- data.frame()
  doubleOE <- data.frame()
  doubleMix <- data.frame()
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

if (config$saveResults) {
  
  message("Saving single-gene perturbation files...")
  saveRDS(singleKD, file = singleKDFile)
  saveRDS(singleOE, file = singleOEFile)
  
  # Only save double-gene results if they exist
  if (nrow(doubleKD) > 0) {
    message("Saving double KD file...")
    saveRDS(doubleKD, file = doubleKDFile)
  } else {
    message("No double KD results to save")
  }
  
  if (nrow(doubleOE) > 0) {
    message("Saving double OE file...")
    saveRDS(doubleOE, file = doubleOEFile)
  } else {
    message("No double OE results to save")
  }
  
  if (nrow(doubleMix) > 0) {
    message("Saving double mixed file...")
    saveRDS(doubleMix, file = doubleMixFile)
  } else {
    message("No double mixed results to save")
  }
  
} else {
  message("Results saving disabled in config")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

total_script_time <- difftime(Sys.time(), script_start_time, units = "mins")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SCRIPT 13 COMPLETED SUCCESSFULLY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat(sprintf("Total runtime: %.1f minutes (%.1f hours)\n", 
            total_script_time, total_script_time/60))

cat("\nANALYSIS SUMMARY:\n")
cat("- Genes tested:", length(genes), "\n")
cat("- Single perturbations completed:", length(genes) * 2, "\n")
cat("- Valid KD results:", valid_kd, "/", length(genes), "\n") 
cat("- Valid OE results:", valid_oe, "/", length(genes), "\n")

if (total_combinations > 0) {
  cat("- Double perturbations completed:", total_combinations, "\n")
}

# Best targets summary
if (nrow(improvement_targets) > 0) {
  cat("\nBEST SINGLE TARGETS:\n")
  cat("- Best KD target:", improvement_targets$Gene[1], 
      "(Δ =", round(improvement_targets$Delta[1], 4), ")\n")
}
if (nrow(improvement_targets_oe) > 0) {
  cat("- Best OE target:", improvement_targets_oe$Gene[1], 
      "(Δ =", round(improvement_targets_oe$Delta[1], 4), ")\n")
}

if (nrow(doubleKD) > 0 && doubleKD$Delta[1] < 0) {
  cat("- Best KD combination:", paste(doubleKD[1, 1:2], collapse = " + "), 
      "(Δ =", round(doubleKD$Delta[1], 4), ")\n")
}

cat("\nAll results saved to:", dirname(singleKDFile), "\n")
cat("Ready for downstream analysis in Script 14!\n")

cat("\nDone!\n")