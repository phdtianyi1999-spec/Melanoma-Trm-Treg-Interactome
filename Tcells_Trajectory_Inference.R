library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
getwd()
setwd("D:/Data_us/data/combined_data/scvi_export_for_r/")

t_cells_clean <- t_cells_clean_final_annotated3
rm(t_cells_clean_final_annotated3)
cd8_obj <- subset(t_cells_clean, subset = major_lineage %in% c("CD8+", "CD8+ (Prolif)"))
cd4_obj <- subset(t_cells_clean, subset = major_lineage == "CD4+")


# Define color palettes
cols_cd4 <- c(
  "CD4+ Naive"                    = "#74C476", 
  "CD4+ Tcm"                      = "#3C5488", 
  "CD4+ T_ISG (STAT1+)"           = "#631879", 
  "CD4+ Treg (Naive)"             = "#00A087", 
  "CD4+ Treg (Intermediate)"      = "#91D1C2", 
  "CD4+ Treg (Activated)"         = "#E64B35", 
  "CD4+ Tfh (CXCL13+)"            = "#B09C85", 
  "CD4+ Th17"                     = "#F39B7F", 
  # Ensure names here exactly match those in the table
  "CD4+ Treg (Th1-Like Transitional)" = "#4DBBD5", 
  "CD4+ Tem"                      = "#BB0021", 
  "CD4+ CTL (NKG7+)"              = "#DC0000",
  "CD4+ Tcm (STAT1+)"             = "#8C6BB1"
)

cols_cd8 <- c(
  # --- Naive & Stem-like ---
  "CD8+ Naive"                  = "#78C679", # Corresponds to Cluster 7
  "CD8+ Tcm (Stem-like)"        = "#41AB5D", # Corresponds to Cluster 1
  
  # --- Exhaustion-related (Cluster 8 & 0) ---
  "CD8+ Tem (GZMK+)"            = "#D9F0A3", # Corresponds to Cluster 8 (Original Pre-Ex color)
  "CD8+ Tex (CXCL13+)"          = "#4A1486", # Corresponds to Cluster 0 (Original Terminal color)
  
  # --- Effector & Transitional (Cluster 4 & 6) ---
  "CD8+ Pre-Trm"                = "#FD8D3C", # Corresponds to Cluster 4 (Original Transitional color)
  "CD8+ Tef (GZMB+)"            = "#F16913", # Corresponds to Cluster 6
  
  # --- Resident (Cluster 2 & 5) ---
  "CD8+ Trm"                    = "#6A51A3", # Corresponds to Cluster 2 (Original Trm Exhausted color)
  # --- Proliferating ---
  "CD8+ Proliferating"          = "#D9D9D9"  # Corresponds to Cluster 3
)
cols_other <- c("Gd-T"          = "#008B45")

# Combine
master_palette <- c(cols_cd4, cols_cd8, cols_other)



# 1. Double-check data preparation (to prevent environment changes)
# Assuming your CD4 object is named cd4_obj
# Extract data to build CDS
data <- GetAssayData(cd4_obj, assay = "RNA", slot = "counts")
cell_metadata <- cd4_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

# ===================================================
# 2. Core Fix: Synchronize detailed annotation information
# ===================================================
# [Key] Synchronize "Detailed_Type" from Seurat to the Monocle object
colData(cds)$Detailed_Type <- cd4_obj@meta.data[rownames(colData(cds)), "Detailed_Type"]

# ===================================================
# 3. Monocle Native Dimensionality Reduction and Trajectory Construction
# ===================================================
# A. Preprocessing (PCA/LSI)
#    num_dim determines the number of principal components for downstream analysis, 50 is a common default
# Core Fix: Fix random seed
# ===================================================
set.seed(123) # Can be any number, as long as it's consistent
cds <- preprocess_cds(cds, num_dim = 50)

# B. Dimensionality Reduction (UMAP) -> [Key Modification] Use coordinates calculated natively by Monocle
#    Stop importing Seurat coordinates to preserve continuous features of the data
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# C. Clustering (Cluster)
#    Monocle needs to cluster (partition) first to draw trajectory lines
cds <- cluster_cells(cds)

# D. Learn Trajectory (Learn Graph)
#    use_partition = TRUE means if cell differences are huge, trajectories will break into disconnected parts
cds <- learn_graph(cds, use_partition = TRUE)

# ===================================================
# 4. Plotting: View cell type distribution on the Monocle trajectory
# ===================================================
# The UMAP shape here will differ from Seurat's, but it better reflects developmental relationships
p1 <- plot_cells(cds, 
                 color_cells_by          = "Detailed_Type", # Use the synchronized column name
                 label_groups_by_cluster = FALSE, 
                 label_leaves            = FALSE, 
                 label_branch_points     = FALSE,
                 graph_label_size        = 1.5,
                 cell_size               = 0.7) +           # Slightly adjust point size
  ggtitle("A. Monocle3 Native Trajectory") +
  theme(legend.position = "right")

print(p1)

# ===================================================
# Monocle 3 Plot Beautification (Match Seurat Style)
# ===================================================

# 1. Redefine your CD4 specific palette (ensure color consistency)
# 2. Plotting: Adjust font size + apply custom colors
p_traj <- plot_cells(cds, 
                     color_cells_by = "Detailed_Type", # Color by subdivided subgroups
                     label_groups_by_cluster = FALSE,  # Turn off Cluster labels
                     label_cell_groups = TRUE,         # Turn on cell type labels
                     label_leaves = FALSE,             # Turn off leaf labels
                     label_branch_points = FALSE,      # Turn off branch point labels
                     group_label_size = 4,             # [Key] Enlarge label font size (default is 2)
                     cell_size = 0.8,                    # Make points slightly larger for better texture
                     trajectory_graph_segment_size = 0.8, # Make the trajectory line slightly thicker
                     graph_label_size = 3) +           # Trajectory node label size
  
  # [Key] Force apply your Nature color palette
  scale_color_manual(values = cols_cd4) +
  
  # Remove cluttered background to keep it clean
  theme_void() +
  ggtitle("Trajectory on CD4 Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right") # Show legend

print(p_traj)
ggsave("Figure_Monocle_Trajectory_Fixed.pdf", width = 10, height = 8)

# Interactively select the starting point (Naive)
cds <- order_cells(cds)

# Plot pseudotime graph (no custom palette needed, uses yellow-purple scale for time)
p_pseudo <- plot_cells(cds, 
                       color_cells_by = "pseudotime", 
                       label_cell_groups = FALSE, 
                       label_leaves = FALSE, 
                       label_branch_points = FALSE,
                       cell_size = 0.8) +
  theme_void() +
  ggtitle("Evolutionary Time (Pseudotime)")

print(p_pseudo)
ggsave("Figure_Monocle_Pseudotime_Fixed.pdf", width = 10, height = 8)

# ===================================================
# 4. Interactively define starting point (Order Cells)
# ===================================================
# After running this line, RStudio will pop up a new window (or in the Plot pane)
# Please click the "circle" or "node" where "CD4+ Naive" cells are most concentrated, then click "Done"
cds <- order_cells(cds)

# ===================================================
# 5. Plotting B: Pseudotime Coloring (Pseudotime)
# ===================================================
# Brighter colors (yellow/bright blue) represent later stages of differentiation
p2 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 graph_label_size = 1.5) +
  ggtitle("B. Pseudotime (Evolutionary Time)")

print(p2)
ggsave("Figure_Monocle_Pseudotime.pdf", width = 8, height = 6)


# ===================================================
# 0. Preparation: Load packages and define CD8 palette
# ===================================================
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

# ===================================================
# CD8 Subgroup Pseudotime Analysis 
# 1. Build Monocle Object (Extract from cd8_obj)
# Assuming cd8_obj is already in your environment
cd8_obj <- subset(t_cells_clean, subset = major_lineage %in% c("CD8+", "CD8+ (Prolif)"))

# Extract Counts matrix (using counts from RNA assay)
data <- GetAssayData(cd8_obj, assay = "RNA", slot = "counts")
cell_metadata <- cd8_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create CDS object
cds_cd8 <- new_cell_data_set(data,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)

# [Key Sync] Synchronize Detailed_Type info to CDS for subsequent plotting
colData(cds_cd8)$Detailed_Type <- cd8_obj@meta.data[rownames(colData(cds_cd8)), "Detailed_Type"]

# ===================================================
# 2. Monocle Native Dimensionality Reduction and Trajectory Construction (without using Seurat coords)
# ===================================================
# A. Preprocessing
set.seed(123) # Set seed to ensure reproducible results
cds_cd8 <- preprocess_cds(cds_cd8, num_dim = 50)

# B. Dimensionality Reduction (Monocle Native UMAP)
# This recalculates coordinates based on data features, usually better suited for pseudotime than Seurat coords
cds_cd8 <- reduce_dimension(cds_cd8, reduction_method = "UMAP")

# C. Clustering (Partition)
cds_cd8 <- cluster_cells(cds_cd8)

# D. Learn Trajectory (Learn Graph)
# use_partition = TRUE allows trajectories to break (e.g., when composed of two completely unrelated cell populations)
cds_cd8 <- learn_graph(cds_cd8, use_partition = TRUE)

# ===================================================
# 3. Plotting A: Cell Type Distribution (Trajectory Plot)
# ===================================================
p_traj_cd8 <- plot_cells(cds_cd8, 
                         color_cells_by = "Detailed_Type", # Color by subdivided subgroups
                         label_groups_by_cluster = FALSE,  
                         label_cell_groups = TRUE,         # Show cell type labels
                         label_leaves = FALSE,             
                         label_branch_points = FALSE,      
                         group_label_size = 4,             # Enlarge font
                         cell_size = 0.8,                  # Point size
                         trajectory_graph_segment_size = 0.8, # Trajectory line thickness
                         graph_label_size = 3) +           
  
  # Apply CD8 specific palette
  scale_color_manual(values = cols_cd8) +
  
  theme_void() +
  ggtitle("Trajectory of CD8+ T Cells") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right")

print(p_traj_cd8)
ggsave("Figure_Monocle_CD8_Trajectory.pdf", width = 10, height = 8)

# ===================================================
# 4. Interactively define starting point (Define Pseudotime)
# ===================================================
# [Operation Guide]
# After running the line below, in the pop-up window please:
# 1. Find the area where "CD8+ Naive" clusters (usually light green points).
# 2. Click the "circle node" extending from this area as the starting point.
# 3. Click "Done".
cds_cd8 <- order_cells(cds_cd8)

# ===================================================
# 5. Plotting B: Pseudotime Coloring (Pseudotime)
# ===================================================
p_pseudo_cd8 <- plot_cells(cds_cd8, 
                           color_cells_by = "pseudotime", 
                           label_cell_groups = FALSE, 
                           label_leaves = FALSE, 
                           label_branch_points = FALSE,
                           cell_size = 0.8) +
  theme_void() +
  ggtitle("CD8+ Evolutionary Time (Pseudotime)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(p_pseudo_cd8)
ggsave("Figure_Monocle_CD8_Pseudotime.pdf", width = 8, height = 6)

# Pseudotime Gene Dynamic Changes Analysis
library(monocle3)
library(ggplot2)

# 1. Set starting point (if not done yet)
# This step is crucial, time must flow from Naive towards Th1-like
# Use helper function to select the node containing Naive as root
# === Define helper function: Find the node closest to the specified cell population ===
get_earliest_principal_node <- function(cds, cell_type_column = "final_cellchat_group", root_type = "CD4+ Naive"){
  
  # 1. Find cell IDs belonging to the "root" type
  cell_ids <- which(colData(cds)[, cell_type_column] == root_type)
  
  # Check if cells were found
  if(length(cell_ids) == 0) {
    stop(paste("Error: Cell named", root_type, "not found in column", cell_type_column, "! Please check spelling."))
  }
  
  # 2. Get the closest nodes to these cells in the Principal Graph
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  # 3. Find which node "captured" the most root cells; that node is the Root Node
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}
# 1. Set starting point (Root)
# Assuming your annotation column is named "final_cellchat_group" and the starting point is named "CD4+ Naive"
# If your starting point is "Treg_Naive", change root_type to "Treg_Naive"
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, 
                                                                    cell_type_column = "Detailed_Type", 
                                                                    root_type = "CD4+ Naive"))

# 2. Check if pseudotime was calculated successfully
# Check if the pseudotime column exists in metadata, or just draw a plot directly
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

# 3. Define key genes (Signaling -> Function)
genes_signaling <- c("STAT1", "IL12RB2", "TNFAIP3", "IRF1")  # Mid-stage: Respond to signals
genes_function  <- c("BATF", "CXCL13", "PRDM1", "TBX21")     # Late-stage: Effector function

genes_to_plot <- c(genes_signaling, genes_function)

# 4. Extract data
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_to_plot, ]

# 5. Plotting (Fitted Curves)
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_to_plot, ]

# [Key Step] Create a simplified grouping column
# The names here must match your real names in Detailed_Type exactly
# Based on your uploaded graph, I guess the names might be these; please confirm with your table(cds$Detailed_Type)
key_types <- c("CD4+ Treg (Naive)", 
               "CD4+ Treg (Th1-Like Transitional)", 
               "CD4+ Treg (Activated)") # Or "Treg_Effector_Th1"

# Create a new column called plot_group
colData(lineage_cds)$plot_group <- as.character(colData(lineage_cds)$Detailed_Type)

# Change anything not in key_types to "Others"
colData(lineage_cds)$plot_group[!colData(lineage_cds)$plot_group %in% key_types] <- "Others"

# Force set factor level order (Naive -> Trans -> Effector -> Others)
# This ensures the legend order is plotted correctly
colData(lineage_cds)$plot_group <- factor(colData(lineage_cds)$plot_group, 
                                          levels = c(key_types, "Others"))

# 3. Plotting (Now 4 colors are enough)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by = "plot_group", # Use our newly created column
                         min_expr = 0.5, 
                         ncol = 2, 
                         cell_size = 1.5, # Slightly enlarge points to see clearly
                         trend_formula = "~ splines::ns(pseudotime, df=3)") +
  
  # Corresponding colors: Naive(Green), Trans(Yellow), Effector(Red), Others(Gray)
  scale_color_manual(values = c("#ABDDA4", "#FDAE61", "#D7191C", "#E0E0E0")) + 
  
  ggtitle("Sequential Activation: Signaling precedes Function")