# ==============================================================================
# 0. Environment Setup and Data Loading
# ==============================================================================
library(nichenetr)
library(Seurat)
library(tidyverse)

# Set the folder path where your .rds files are stored
data_dir <- "nichenet_db/" 

# Assuming you are working with human data (if mouse, replace 'human' with 'mouse' in the filenames below, and append the '_mouse' suffix)
organism <- "human" 

print("Loading NicheNet knowledge base (files are large, please be patient)...")

# 1. Core prediction model (Ligand -> Target)
ligand_target_matrix <- readRDS(paste0(data_dir, "ligand_target_matrix_nsga2r_final.rds"))

# 2. Ligand-receptor network (Who binds to whom)
lr_network <- readRDS(paste0(data_dir, "lr_network_human_21122021.rds"))

# 3. [Advanced] Ligand-Transcription Factor matrix (Ligand -> TF, used to link with SCENIC)
ligand_tf_matrix <- readRDS(paste0(data_dir, "ligand_tf_matrix_nsga2r_final.rds"))

# 4. Weighted networks (used for visualizing signaling pathways)
weighted_networks <- readRDS(paste0(data_dir, "weighted_networks_nsga2r_final.rds"))

print("✅ All databases loaded! Ready to start analysis...")

# ==============================================================================
# 1. Define Sender and Receiver
# ==============================================================================

# --- Please modify here ---
# Define your cell type names (must match exactly with those in your Seurat object)
receiver_type <- "T_cells"       # Receiver: Your T cells
sender_type <- c("Macrophages")  # Sender: Macrophages

# 1.1 Extract expressed genes in the receiver
expressed_genes_receiver <- get_expressed_genes(receiver_type, seurat_obj, pct = 0.10)

# 1.2 Extract expressed genes in the sender
list_expressed_genes_sender <- sender_type %>% unique() %>% lapply(get_expressed_genes, seurat_obj, pct = 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

# 1.3 Define the Gene Set of interest
# Assuming these are the Top genes found in your NMF glycolysis module
# If you haven't defined them yet, here are a few manual examples:
genes_of_interest <- c("HK2", "LDHA", "PKM", "HIF1A", "GAPDH", "ENO1", "SLC2A1", "PFKP")

# Filter: Must be genes that are actually expressed in T cells
genes_of_interest <- genes_of_interest %>% intersect(expressed_genes_receiver)
background_expressed_genes <- expressed_genes_receiver %>% dplyr::setdiff(genes_of_interest)

# ==============================================================================
# 2. Filter Potential Ligands
# ==============================================================================

# Logic: Sender has the ligand, Receiver has the receptor, and they can bind according to the database
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()

print(paste0("🔍 Found ", length(potential_ligands), " potential interacting ligands."))

# ==============================================================================
# 3. Core calculation: Predict which ligands drive the glycolysis genes?
# ==============================================================================

print("Calculating ligand activities (Pearson Correlation)...")
ligand_activities <- predict_ligand_activities(
  genes_of_interest = genes_of_interest,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Select the top 20 ranked ligands
best_upstream_ligands <- ligand_activities %>% 
  arrange(desc(pearson)) %>% 
  top_n(20, pearson) %>% 
  pull(test_ligand) %>% unique()

print("🏆 Most likely driving ligands (Top Ligands):")
print(best_upstream_ligands)

# ==============================================================================
# 4. Visualization I: Ligand-Target Genes heatmap
# ==============================================================================
# This plot shows: Does Macrophage IL6 (ligand) regulate T cell HK2 (target gene)?

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_matrix = ligand_target_matrix,
  targets_of_interest = genes_of_interest,
  top_ligands = best_upstream_ligands
)

# Plotting
p_ligand_target <- make_heatmap_ggplot(
  active_ligand_target_links, 
  y_name = "Prioritized Ligands (Macrophage)", 
  x_name = "Target Genes (Glycolysis)", 
  color = "darkorange", 
  legend_position = "top", 
  x_axis_position = "top",
  legend_title = "Regulatory Potential"
)
print(p_ligand_target)

# ==============================================================================
# 5. [Advanced] Visualization II: Ligand-TF heatmap (Ligand -> TF)
# ==============================================================================
# This step uses the specifically downloaded ligand_tf_matrix!
# It helps link CellChat and SCENIC together.

# Assuming your SCENIC analysis found high activity for these TFs in T cells:
scenic_tfs <- c("HIF1A", "STAT3", "MYC", "NFKB1", "JUN") 

# Calculate the regulatory potential of ligands on these TFs
active_ligand_tf_links <- prepare_ligand_target_visualization(
  ligand_target_matrix = ligand_tf_matrix, # Note: Switched to the TF matrix here
  targets_of_interest = scenic_tfs,
  top_ligands = best_upstream_ligands
)

# Plotting
p_ligand_tf <- make_heatmap_ggplot(
  active_ligand_tf_links, 
  y_name = "Prioritized Ligands (Macrophage)", 
  x_name = "Transcription Factors (SCENIC)", 
  color = "purple", 
  legend_position = "top", 
  x_axis_position = "top",
  legend_title = "Signaling Potential"
)
print(p_ligand_tf)

# ==============================================================================
# 6. Visualization III: Ligand-Receptor expression dot plot (Validation)
# ==============================================================================
# See if these inferred ligands and receptors are actually highly expressed?

# Find the corresponding receptors for the Top ligands
lr_network_top <- lr_network %>% 
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from, to)
best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

# Plotting
p_dotplot <- DotPlot(seurat_obj, 
                     features = unique(c(best_upstream_ligands, best_upstream_receptors)), 
                     cols = "RdYlBu") + 
  RotatedAxis() +
  ggtitle("Expression of Top Ligands & Receptors")

print(p_dotplot)