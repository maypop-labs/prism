# ======================================================================
# 13 - Graphs
# Purpose: Load saved data objects and generate graphs to illustrate
#          key steps of the project without revealing the actual gene
#          names.
# ======================================================================

# --- Load Required Libraries ---
library(Seurat)
library(monocle3)
library(ggplot2)
library(igraph)
library(dplyr)
library(tidyr)
library(readr)

# --- 1. PCA/UMAP Plot from Seurat Object ---
# Load your merged Seurat object (adjust file path as needed)
seurat_file <- "E:/datasets/omics/skin/results/rds/merged_seurat.rds"  # Adjust path if necessary
seurat_obj <- readRDS(seurat_file)

# Generate a UMAP plot colored by donor age (or any other metadata)
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "age") +
  ggtitle("UMAP Plot of All Cells (by Age)") +
  theme_minimal()

# Save UMAP plot
ggsave("D:/projects/Gene Regulatory Network Perterbation in Aging/R Plots/umap_plot.png", plot = umap_plot, width = 8, height = 6)


# --- 2. Pseudotime Trajectory from Monocle3 Object ---
# Load the Monocle3 object for keratinocytes trajectory ("Y_447" smoothed)
# (Assumes you have a function called load_monocle_objects; adjust as needed.)
monocle_dir <- "E:/datasets/omics/skin/results/monocle3/monocle3_Keratinocytes_Y_447_smoothed"  # Update path if necessary
cds <- load_monocle_objects(directory_path = monocle_dir)  # Custom helper function

# Plot the pseudotime trajectory with cells colored by pseudotime value
pt_plot <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE
) + 
  ggtitle("Pseudotime Trajectory") +
  theme_minimal()

# Save pseudotime plot
ggsave("D:/projects/Gene Regulatory Network Perterbation in Aging/R Plots/pseudotime_plot.png", plot = pt_plot, width = 8, height = 6)


# --- 3. Gene Regulatory Network (GRN) Visualization ---
# Load the GRN edges RDS file
grn_edges_file <- "E:/datasets/omics/skin/results/rds/Keratinocytes_Y_447_GRN_Part_02_edges.rds"
grn_edges <- readRDS(grn_edges_file)

# Expect grn_edges to have columns "TF" and "Target". Create a mapping to generic names.
all_genes <- unique(c(as.character(grn_edges$TF), as.character(grn_edges$Target)))
generic_names <- paste0("Gene", seq_along(all_genes))
gene_mapping <- setNames(generic_names, all_genes)

# Replace gene names with generic names in the edge table
grn_edges$TF_gen <- gene_mapping[as.character(grn_edges$TF)]
grn_edges$Target_gen <- gene_mapping[as.character(grn_edges$Target)]

# Build an igraph network from the pruned edge table
grn_graph <- graph_from_data_frame(grn_edges[, c("TF_gen", "Target_gen")], directed = TRUE)

# Plot the network. Adjust vertex size and arrow size as desired.
png("grn_network.png", width = 800, height = 600)
plot(grn_graph,
     vertex.label = V(grn_graph)$name,
     vertex.size = 5,
     vertex.label.cex = 0.8,
     edge.arrow.size = 0.5,
     main = "Gene Regulatory Network (Generic Gene Names)")
dev.off()


# --- 4. Attractor Landscape Plot ---
# Load attractor data frame containing AgeScore and BasinSize (anonymized)
attractor_df_file <- "E:/datasets/omics/skin/results/rds/Keratinocytes_Y_447_attractor_df.rds"
attractor_df <- readRDS(attractor_df_file)
# Expect attractor_df to have columns such as Attractor, AgeScore, BasinSize.
# Relabel attractors generically (e.g., Attractor1, Attractor2, etc.)
attractor_df$Attractor <- paste0("Attractor", attractor_df$Attractor)

landscape_plot <- ggplot(attractor_df, aes(x = AgeScore, y = BasinSize)) +
  geom_point(color = "blue", size = 3) +
  #geom_text(aes(label = Attractor), vjust = -0.5) +
  labs(title = "Attractor Landscape", x = "Age Score (0 = Young, 1 = Old)", y = "Basin Size Fraction") +
  theme_minimal()

# Save attractor landscape plot
ggsave("D:/projects/Gene Regulatory Network Perterbation in Aging/R Plots/attractor_landscape.png", plot = landscape_plot, width = 8, height = 6)


# --- 5. Perturbation Analysis Results ---
# Load final single-gene perturbation results (two files for different modes, e.g., knockdown and overexpression)
perturb_file0 <- "E:/datasets/omics/skin/results/rds/Keratinocytes_Y_447_final_single_targets_0.rds"
perturb_file1 <- "E:/datasets/omics/skin/results/rds/Keratinocytes_Y_447_final_single_targets_1.rds"
perturb_data0 <- readRDS(perturb_file0)
perturb_data1 <- readRDS(perturb_file1)

# Combine and label the mode (here we assume a column "Gene" and "AgingScore" exist)
perturb_data0 <- perturb_data0 %>% mutate(Mode = "Knockdown")
perturb_data1 <- perturb_data1 %>% mutate(Mode = "Overexpression")
perturb_df <- bind_rows(perturb_data0, perturb_data1)

# Map actual gene names to generic names
unique_perturb_genes <- unique(as.character(perturb_df$Gene))
perturb_generic <- paste0("Gene", seq_along(unique_perturb_genes))
perturb_mapping <- setNames(perturb_generic, unique_perturb_genes)
perturb_df$Gene_gen <- perturb_mapping[as.character(perturb_df$Gene)]

# Create a dot plot of Aging Score per gene for each perturbation mode
perturb_plot <- ggplot(perturb_df, aes(x = reorder(Gene_gen, AgingScore), y = AgingScore, color = Mode)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(title = "Perturbation Analysis: Aging Scores (Generic Genes)",
       x = "Gene (Generic)",
       y = "Aging Score") +
  theme_minimal()

# Save perturbation plot
ggsave("D:/projects/Gene Regulatory Network Perterbation in Aging/R Plots/perturbation_results.png", plot = perturb_plot, width = 8, height = 6)


# --- 6. (Optional) PCA Plot from Seurat Object ---
# In case you want a PCA plot as well
pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "age") +
  ggtitle("PCA Plot (By Age)") +
  theme_minimal()
ggsave("D:/projects/Gene Regulatory Network Perterbation in Aging/R Plots/pca_plot.png", plot = pca_plot, width = 8, height = 6)

# ======================================================================
# End of Script
# ======================================================================
