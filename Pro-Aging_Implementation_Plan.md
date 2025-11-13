# Pro-Aging (Anti-Hit) Implementation Plan for Script 08

## Overview
Add tracking and analysis of perturbations that **increase** the attractor landscape aging score (Delta > 0.01). These "pro-aging" or "anti-hit" perturbations represent genes whose perturbation accelerates aging and should be avoided therapeutically.

## Changes Required

### 1. Part 3 - Add Pro-Aging Counters

After line 250 (where we count improvements), add pro-aging counters:

```r
# Count pro-aging effects (anti-hits)
kdProAging <- sum(singleKD$Delta > 0.01, na.rm = TRUE)
oeProAging <- sum(singleOE$Delta > 0.01, na.rm = TRUE)

message("  - KD pro-aging effects: ", kdProAging)
message("  - OE pro-aging effects: ", oeProAging)
```

Update `analysisSummary` to include pro-aging counts:

```r
analysisSummary <- data.frame(
  InitialAgingScore = initialScore,
  GenesAnalyzed = nGenes,
  ValidKDResults = validKD,
  ValidOEResults = validOE,
  KDImprovements = kdImprovements,
  OEImprovements = oeImprovements,
  KDProAging = kdProAging,
  OEProAging = oeProAging,
  stringsAsFactors = FALSE
)
```

### 2. Part 4 - Extract and Rank Pro-Aging Hits

After the improvements section (around line 320), add:

```r
# Extract pro-aging effects (anti-hits)
proAgingHits <- successful[successful$Delta > 0.01, ]

message("Pro-aging effects identified:")
message("  - Total pro-aging hits: ", nrow(proAgingHits))

# Initialize empty object
proAgingPivotTable <- data.frame()

if (nrow(proAgingHits) > 0) {
  # Sort by worsening (most positive delta first)
  proAgingHits <- proAgingHits[order(proAgingHits$Delta, decreasing = TRUE), ]
  
  # Add ranking
  proAgingHits$Rank <- 1:nrow(proAgingHits)
  proAgingHits$PercentWorsening <- abs(proAgingHits$Delta) * 100 / abs(proAgingHits$Delta[1])
  
  message("\nTop 5 pro-aging effects to AVOID:")
  topProAgingResults <- head(proAgingHits[, c("Gene", "Mode", "Delta", "Rank")], 5)
  print(topProAgingResults)
  
  # Create pivot table for pro-aging hits
  proAgingData <- proAgingHits[, c("Gene", "Mode", "Delta", "AgingScore")]
  
  kdProAgingData <- proAgingData[proAgingData$Mode == "Knockdown", c("Gene", "Delta", "AgingScore")]
  oeProAgingData <- proAgingData[proAgingData$Mode == "Overexpression", c("Gene", "Delta", "AgingScore")]
  
  colnames(kdProAgingData) <- c("Gene", "KD_Delta", "KD_AgingScore")
  colnames(oeProAgingData) <- c("Gene", "OE_Delta", "OE_AgingScore")
  
  proAgingPivotTable <- merge(kdProAgingData, oeProAgingData, by = "Gene", all = TRUE)
  
  # Calculate worst result per gene
  proAgingPivotTable$WorstDelta <- pmax(proAgingPivotTable$KD_Delta, proAgingPivotTable$OE_Delta, na.rm = TRUE)
  proAgingPivotTable$WorstMode <- ifelse(is.na(proAgingPivotTable$KD_Delta), "OE",
                                  ifelse(is.na(proAgingPivotTable$OE_Delta), "KD",
                                  ifelse(proAgingPivotTable$KD_Delta > proAgingPivotTable$OE_Delta, "KD", "OE")))
  
  # Sort by worst result
  proAgingPivotTable <- proAgingPivotTable[order(proAgingPivotTable$WorstDelta, decreasing = TRUE), ]
  
  message("Part 4 completed: ", nrow(proAgingPivotTable), " pro-aging genes identified")
}
```

### 3. Part 5 - Track Pro-Aging Double Perturbations

After counting double improvements (line 455), add:

```r
# Count pro-aging double effects
doubleKDProAging <- sum(doubleKD$Delta > 0.01, na.rm = TRUE)
doubleOEProAging <- sum(doubleOE$Delta > 0.01, na.rm = TRUE)
doubleMixProAging <- sum(doubleMix$Delta > 0.01, na.rm = TRUE)

message("  - Double KD pro-aging: ", doubleKDProAging)
message("  - Double OE pro-aging: ", doubleOEProAging)
message("  - Mixed pro-aging: ", doubleMixProAging)
```

Find worst double perturbation at the end of Part 5:

```r
# Find worst double perturbation overall
worstDoubleScore <- NA_real_
worstDoubleDelta <- NA_real_
worstDoubleType <- "None"
worstDoubleGenes <- "None"

if (nrow(doubleKD) > 0) {
  # Sort by delta descending to get worst at top
  doubleKDSorted <- doubleKD[order(doubleKD$Delta, decreasing = TRUE), ]
  if (!is.na(doubleKDSorted$Delta[1]) && doubleKDSorted$Delta[1] > 0.01) {
    worstDoubleScore <- doubleKDSorted$AgingScore[1]
    worstDoubleDelta <- doubleKDSorted$Delta[1]
    worstDoubleType <- "Double KD"
    worstDoubleGenes <- paste(doubleKDSorted$Gene1[1], "+", doubleKDSorted$Gene2[1])
  }
}

if (nrow(doubleOE) > 0) {
  doubleOESorted <- doubleOE[order(doubleOE$Delta, decreasing = TRUE), ]
  if (!is.na(doubleOESorted$Delta[1]) && doubleOESorted$Delta[1] > 0.01 &&
      (is.na(worstDoubleDelta) || doubleOESorted$Delta[1] > worstDoubleDelta)) {
    worstDoubleScore <- doubleOESorted$AgingScore[1]
    worstDoubleDelta <- doubleOESorted$Delta[1]
    worstDoubleType <- "Double OE"
    worstDoubleGenes <- paste(doubleOESorted$Gene1[1], "+", doubleOESorted$Gene2[1])
  }
}

if (nrow(doubleMix) > 0) {
  doubleMixSorted <- doubleMix[order(doubleMix$Delta, decreasing = TRUE), ]
  if (!is.na(doubleMixSorted$Delta[1]) && doubleMixSorted$Delta[1] > 0.01 &&
      (is.na(worstDoubleDelta) || doubleMixSorted$Delta[1] > worstDoubleDelta)) {
    worstDoubleScore <- doubleMixSorted$AgingScore[1]
    worstDoubleDelta <- doubleMixSorted$Delta[1]
    worstDoubleType <- "Mixed KD/OE"
    worstDoubleGenes <- paste(doubleMixSorted$KD_Gene[1], "(KD) +", doubleMixSorted$OE_Gene[1], "(OE)")
  }
}

if (!is.na(worstDoubleDelta)) {
  message("\nWorst double perturbation (AVOID):")
  message("  Type: ", worstDoubleType)
  message("  Genes: ", worstDoubleGenes)
  message("  Worsening: +", round(worstDoubleDelta, 4))
}
```

### 4. Part 6 - Add Switch-Gene Pro-Aging Tracking

After switch-gene improvements counting (line 597), add:

```r
switchKDProAging <- sum(switchSingleKD$Delta > 0.01, na.rm = TRUE)
switchOEProAging <- sum(switchSingleOE$Delta > 0.01, na.rm = TRUE)

message("  - KD pro-aging effects: ", switchKDProAging)
message("  - OE pro-aging effects: ", switchOEProAging)
```

After switch double counting (line 716), add:

```r
switchDoubleKDProAging <- sum(switchDoubleKD$Delta > 0.01, na.rm = TRUE)
switchDoubleOEProAging <- sum(switchDoubleOE$Delta > 0.01, na.rm = TRUE)
switchDoubleMixProAging <- sum(switchDoubleMix$Delta > 0.01, na.rm = TRUE)

message("  - Double KD pro-aging: ", switchDoubleKDProAging)
message("  - Double OE pro-aging: ", switchDoubleOEProAging)
message("  - Mixed pro-aging: ", switchDoubleMixProAging)
```

Extract switch pro-aging hits:

```r
switchProAgingHits <- switchSuccessful[switchSuccessful$Delta > 0.01, ]

if (nrow(switchProAgingHits) > 0) {
  switchProAgingHits <- switchProAgingHits[order(switchProAgingHits$Delta, decreasing = TRUE), ]
  message("\nTop 5 pro-aging switch-gene effects to AVOID:")
  topSwitchProAging <- head(switchProAgingHits[, c("Gene", "Mode", "Delta")], 5)
  print(topSwitchProAging)
}
```

### 5. Save Section - Add Pro-Aging Outputs

In the save section, add after Part 4 outputs:

```r
# Pro-aging results (full network)
if (nrow(proAgingHits) > 0) {
  saveObject(proAgingHits, ptPaths$proAgingSummaryTsv, config, "pro-aging summary (TSV)")
  
  # Separate single and double pro-aging
  proAgingSingleKD <- singleKD[singleKD$Success & !is.na(singleKD$Delta) & singleKD$Delta > 0.01, ]
  proAgingSingleOE <- singleOE[singleOE$Success & !is.na(singleOE$Delta) & singleOE$Delta > 0.01, ]
  
  if (nrow(proAgingSingleKD) > 0) {
    saveObject(proAgingSingleKD, ptPaths$proAgingSingleKD, config, "pro-aging single KD")
  }
  if (nrow(proAgingSingleOE) > 0) {
    saveObject(proAgingSingleOE, ptPaths$proAgingSingleOE, config, "pro-aging single OE")
  }
}

# Pro-aging double perturbations
proAgingDoubleKD <- doubleKD[!is.na(doubleKD$Delta) & doubleKD$Delta > 0.01, ]
proAgingDoubleOE <- doubleOE[!is.na(doubleOE$Delta) & doubleOE$Delta > 0.01, ]
proAgingDoubleMix <- doubleMix[!is.na(doubleMix$Delta) & doubleMix$Delta > 0.01, ]

if (nrow(proAgingDoubleKD) > 0) {
  saveObject(proAgingDoubleKD, ptPaths$proAgingDoubleKD, config, "pro-aging double KD")
}
if (nrow(proAgingDoubleOE) > 0) {
  saveObject(proAgingDoubleOE, ptPaths$proAgingDoubleOE, config, "pro-aging double OE")
}
if (nrow(proAgingDoubleMix) > 0) {
  saveObject(proAgingDoubleMix, ptPaths$proAgingDoubleMix, config, "pro-aging double Mix")
}
```

Add for switch-gene pro-aging (after Part 6 outputs):

```r
# Pro-aging switch-gene results
if (nrow(switchProAgingHits) > 0) {
  # Extract switch pro-aging singles
  proAgingSwitchSingleKD <- switchSingleKD[switchSingleKD$Success & !is.na(switchSingleKD$Delta) & switchSingleKD$Delta > 0.01, ]
  proAgingSwitchSingleOE <- switchSingleOE[switchSingleOE$Success & !is.na(switchSingleOE$Delta) & switchSingleOE$Delta > 0.01, ]
  
  if (nrow(proAgingSwitchSingleKD) > 0) {
    saveObject(proAgingSwitchSingleKD, ptPaths$proAgingSwitchSingleKD, config, "pro-aging switch single KD")
  }
  if (nrow(proAgingSwitchSingleOE) > 0) {
    saveObject(proAgingSwitchSingleOE, ptPaths$proAgingSwitchSingleOE, config, "pro-aging switch single OE")
  }
  
  # Extract switch pro-aging doubles
  proAgingSwitchDoubleKD <- switchDoubleKD[!is.na(switchDoubleKD$Delta) & switchDoubleKD$Delta > 0.01, ]
  proAgingSwitchDoubleOE <- switchDoubleOE[!is.na(switchDoubleOE$Delta) & switchDoubleOE$Delta > 0.01, ]
  proAgingSwitchDoubleMix <- switchDoubleMix[!is.na(switchDoubleMix$Delta) & switchDoubleMix$Delta > 0.01, ]
  
  if (nrow(proAgingSwitchDoubleKD) > 0) {
    saveObject(proAgingSwitchDoubleKD, ptPaths$proAgingSwitchDoubleKD, config, "pro-aging switch double KD")
  }
  if (nrow(proAgingSwitchDoubleOE) > 0) {
    saveObject(proAgingSwitchDoubleOE, ptPaths$proAgingSwitchDoubleOE, config, "pro-aging switch double OE")
  }
  if (nrow(proAgingSwitchDoubleMix) > 0) {
    saveObject(proAgingSwitchDoubleMix, ptPaths$proAgingSwitchDoubleMix, config, "pro-aging switch double Mix")
  }
  
  # Generate comprehensive switch pro-aging summary
  message("Generating switch-gene pro-aging summary...")
  
  allSwitchProAgingResults <- data.frame(
    Rank = integer(),
    Target = character(),
    Type = character(),
    AgingScore = numeric(),
    Delta = numeric(),
    SwitchDirection = character(),
    Pseudotime = numeric(),
    PseudoR2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add formatted results similar to improvement section
  # ... (similar structure to allSwitchTargetsResults but for pro-aging)
  
  if (nrow(allSwitchProAgingResults) > 0) {
    allSwitchProAgingResults <- allSwitchProAgingResults[order(allSwitchProAgingResults$Delta, decreasing = TRUE), ]
    allSwitchProAgingResults$Rank <- 1:nrow(allSwitchProAgingResults)
    
    write.table(allSwitchProAgingResults, ptPaths$proAgingSwitchSummaryTsv, sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("Switch-gene pro-aging summary saved: ", nrow(allSwitchProAgingResults), " total")
  }
}
```

### 6. Final Summary - Add Pro-Aging Statistics

At the end (around line 1050), add after the improvements summary:

```r
# Pro-aging summary
if (nrow(proAgingHits) > 0) {
  worstKD <- max(singleKD$Delta, na.rm = TRUE)
  worstOE <- max(singleOE$Delta, na.rm = TRUE)
  worstSingle <- max(worstKD, worstOE)
  
  message("\nPro-aging effects found (AVOID THESE): ", nrow(proAgingHits))
  message("Worst single pro-aging effect: +", round(worstSingle, 4), " (", round(worstSingle * 100, 2), "% increase)")
  
  if (nrow(proAgingHits) > 0) {
    message("\nWorst single pro-aging target (AVOID): ", proAgingHits$Gene[1], " (", proAgingHits$Mode[1], ")")
    message("  Worsening: +", round(proAgingHits$Delta[1], 4))
  }
}
```

## File Naming Convention

**Pro-aging files** use the `_proaging_` prefix to clearly distinguish them from anti-aging targets:
- `*_proaging_single_KD.rds`
- `*_proaging_single_OE.rds`  
- `*_proaging_double_KD_KD.rds`
- `*_proaging_double_OE_OE.rds`
- `*_proaging_double_KD_OE.rds`
- `*_proaging_summary.tsv`
- `*_proaging_switch_single_KD.rds`
- etc.

## Output Structure

Each pro-aging file will contain the same columns as the improvement files but sorted by **worst (most positive) Delta first**:
- Gene/Target
- Mode/Type
- AgingScore
- Delta (positive values)
- Rank
- PercentWorsening (analogous to PercentImprovement)

This allows researchers to:
1. Identify which perturbations accelerate aging
2. Avoid these in therapeutic design
3. Understand which genes are critical for maintaining youthful states
