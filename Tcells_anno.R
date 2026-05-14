# First load data placed in D drive data_us
library(Seurat)
library(ggplot2)
library(dplyr)
getwd()
setwd("D:/Data_us/data/res/")
# ===================================================
# 1. Define sample group lists

# pPR group (Partial Response)
pPR_samples <- c("FZD-1", "FZD-2", "ZYC-1", "ZYC-2", "XGS-1", "XGS-2", "LYZ-1", "LYZ-2")
# pNR group (No Response)
pNR_samples <- c("QRJ-2", "QRK-1", "CXF-1", "CXF-2", "YZY-1", "YZY-2", "LYH-1", "LYH-2")
# PD group (Progressive Disease)
PD_samples <- c("LSX-1", "ZWX-1", "YYX-1")

# 2. Add Reactivity column

# Set to "Unknown" by default first (covers T series)
sce.final$Reactivity <- "Unknown"

# Assign values sequentially
sce.final$Reactivity[sce.final$orig.ident %in% pPR_samples] <- "pPR"
sce.final$Reactivity[sce.final$orig.ident %in% pNR_samples] <- "pNR"
sce.final$Reactivity[sce.final$orig.ident %in% PD_samples] <- "PD"

# 3. Add Timepoint (Pre/Post treatment) column

# Logic: suffix "-1" is Pre, "-2" is Post
# T series has no suffix -1/-2, keep as NA or Unknown
sce.final$Timepoint <- NA

# Use regex matching: ending with "-1" set to Pre
sce.final$Timepoint[grep("-1$", sce.final$orig.ident)] <- "Pre"

# Ending with "-2" set to Post
sce.final$Timepoint[grep("-2$", sce.final$orig.ident)] <- "Post"

# T series keep as NA, if you confirm T series is also pre-treatment, you can set it manually
# sce.final$Timepoint[grep("^T", sce.final$orig.ident)] <- "Unknown" 

# 4. Add Group column (combined grouping) - Recommended for differential analysis
# Format example: pPR_Pre, pNR_Post
# For T series or incomplete info, this column will become NA_Unknown or similar, won't affect analysis
sce.final$Group <- paste(sce.final$Reactivity, sce.final$Timepoint, sep = "_")

# Clean up T series Group names (optional)
sce.final$Group[grep("Unknown", sce.final$Group)] <- "Unknown"

# 5. Check results

message("✅ Metadata addition complete! Check the grouping statistics:")
print(table(sce.final$Reactivity, sce.final$Timepoint))

message("\nCombined grouping (Group) statistics:")
print(table(sce.final$Group))

# 6. Save the object with complete clinical info
saveRDS(sce.final, "sce_final_clinical.rds")


# ===================================================
# Keep cohort samples with treatment info
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Assuming your previous object is called sce.final
# 1. Remove Unknown (T series) samples, keep only those with clinical info
sce.clinical <- subset(sce.final, subset = Reactivity != "Unknown")

# 2. Add Response (binary classification) column
# Logic: pPR -> Responder (R); pNR/PD -> Non-Responder (NR)
sce.clinical$Response <- ifelse(sce.clinical$Reactivity == "pPR", "Responder", "Non-Responder")

# Check if grouping is correct
table(sce.clinical$Reactivity, sce.clinical$Response)
# Should show: pPR all in Responder, pNR and PD all in Non-Responder

# 3. Extract T cells
# Based on your previous annotation, extract clusters containing the word "T Cell"
# Or directly specify idents: c("CD8+ T Cells", "T Cells", "Tregs", "CD4+ T Cells"...)
# Assuming your cell_type column contains "T Cell" or "Tregs"
t_cells <- subset(sce.clinical, subset = cell_type %in% c("T Cells"))

# 4. T cell Re-clustering
# Because a subset was extracted, it's recommended to rerun dimensionality reduction to better separate T cell subclusters
DefaultAssay(t_cells) <- "RNA"
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells)

# The resolution here can be adjusted slightly higher (e.g., 0.6 or 0.8) because T cell subclusters are very fine
t_cells <- FindNeighbors(t_cells, dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.8)
t_cells <- RunUMAP(t_cells, dims = 1:20)

# Plot to see the current T cell subclusters (you will need to re-annotate what these 0,1,2... T cells are)
DimPlot(t_cells, reduction = "umap", label = TRUE) + ggtitle("Re-clustered T Cells")

# ===================================================
# T cell annotation
# Major lineage classification
# 1. Define lineage differentiation Markers
# CD3D/E: Proves it's a T cell
# CD4: CD4 T cell
# CD8A/CD8B: CD8 T cell
# TRDC: Gamma-delta T cell (prevent contamination)
# MKI67: Proliferating cell (sometimes CD4/8 in proliferating cells express weakly, look specifically)

lineage_markers <- c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRDC", "MKI67","SLC4A10", "KLRB1")

# 2. Plot DotPlot (to see expression proportions), take a quick look to get an idea, then custom cluster
DotPlot(t_cells, features = lineage_markers, group.by = "seurat_clusters") + 
  RotatedAxis()

# 1. Define final major lineage annotation mapping (Start custom clustering, based on previous DotPlot)
t_lineage_ids <- c(
  # --- CD8 Family ---
  "0"  = "CD8+",
  "2"  = "CD8+",
  "3"  = "CD8+",         
  "6"  = "CD8+",         
  "7"  = "CD8+ (Prolif)",
  "8"  = "CD8+",
  "10" = "CD8+",
  "12" = "CD8+",         
  "18" = "CD8+",         
  
  # --- CD4 Family ---
  "1"  = "CD4+",
  "5"  = "CD4+",
  "9"  = "CD4+",
  "13" = "CD4+",         
  "15" = "CD4+",
  "16" = "CD4+",         
  "17" = "CD4+",
  "19" = "CD4+",
  "20" = "CD4+",
  
  # --- Special/Other ---
  "14" = "Gd-T",         # Gamma-Delta
  "11" = "Doublet",      # Double positive
  "4"  = "DN"            # Double negative
)

# 2. Rename identities (If your Idents are already names like CD4+, this step is just overwrite confirmation)
t_cells <- RenameIdents(t_cells, t_lineage_ids)

# 3. [Key Step] Save current Identities into the "major_lineage" column
# This step is to resolve errors
t_cells$major_lineage <- Idents(t_cells)

# 4. Check if successful (Should see a column named major_lineage)
print(colnames(t_cells@meta.data))

# Plot DotPlot for T cell major lineage markers
library(Seurat)
library(ggplot2)

# Define gene list
pub_markers <- c("CD3D", "CD3E",          # T Cell Core
                 "CD4",                   # CD4 Lineage
                 "CD8A", "CD8B",          # CD8 Lineage
                 "TRDC",                  # Gamma-Delta
                 "MKI67", "TOP2A",        # Proliferation
                 "SLC4A10", "KLRB1",      # MAIT / Th17 check
                 "GNLY", "NKG7")          # Cytotoxicity

# Plot DotPlot
p_dot <- DotPlot(t_cells, 
                 features = pub_markers, 
                 group.by = "major_lineage", # Definitely found this time
                 cols = c("#E6E6E6", "#D51F26"), 
                 dot.scale = 8) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  labs(x = "Marker Genes", y = "Major Lineage", title = "Lineage Definition Markers")

print(p_dot)

# Save
ggsave("Figure_T_Cell_Lineage_DotPlot.pdf", width = 8, height = 5)


# Plot UMAP for T cell major lineage annotation
# 1. Ensure current default identity is "major_lineage"
# (If you ran t_cells$major_lineage <- Idents(t_cells) previously, this should be fine)
Idents(t_cells) <- "major_lineage"

# 2. Plot UMAP
# Define precisely matched Morandi colors
# This way, no matter how cell order changes in your plot, CD8 is always blue, CD4 is always pink
morandi_specific <- c(
  "CD8+"          = "#8F9FBC",  # Haze blue (corresponds to cytotoxic/cool color)
  "CD8+ (Prolif)" = "#7A8B9D",  # Darker haze blue
  "CD4+"          = "#D6B0B1",  # Dusty pink (corresponds to helper/warm color)
  "Gd-T"          = "#88A096",  # Grey-green (special T cells)
  "DN"            = "#D4C4B7",  # Milk coffee color (neutral)
  "Doublet"       = "#9E8E8E"   # Grey (to be removed)
)

# 2. Plotting
p_umap_custom <- DimPlot(t_cells, 
                         reduction = "umap", 
                         label = TRUE, 
                         label.size = 5, 
                         repel = TRUE, 
                         pt.size = 0.5, 
                         cols = morandi_specific) + # Seurat automatically matches names
  ggtitle("T Cell Major Lineages") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(p_umap_custom)

# Save
ggsave("Figure_T_Cell_UMAP_Morandi_Custom.pdf", width = 8, height = 7)




# ===================================================
# 1. Data cleaning and rerun UMAP
message("Removing Doublet and DN cells...")
# invert = TRUE means "invert selection", i.e., keep cells other than these
t_cells_clean <- subset(t_cells, idents = c("Doublet", "DN"), invert = TRUE)

# Rerun standard pipeline (Since cell count changed, variance structure changed, recommend rerunning PCA and UMAP)
DefaultAssay(t_cells_clean) <- "RNA"
t_cells_clean <- NormalizeData(t_cells_clean)
t_cells_clean <- FindVariableFeatures(t_cells_clean, nfeatures = 2000)
t_cells_clean <- ScaleData(t_cells_clean)
t_cells_clean <- RunPCA(t_cells_clean, verbose = FALSE)

# Rerun UMAP
t_cells_clean <- RunUMAP(t_cells_clean, dims = 1:20)

# Define Morandi color palette (maintain consistency)
morandi_specific <- c(
  "CD8+"          = "#8F9FBC",  # Haze blue
  "CD8+ (Prolif)" = "#7A8B9D",  # Dark haze blue
  "CD4+"          = "#D6B0B1",  # Dusty pink
  "Gd-T"          = "#88A096"   # Grey-green
)

# Plot cleaned UMAP
p_umap_clean <- DimPlot(t_cells_clean, 
                        reduction = "umap", 
                        cols = morandi_specific,
                        label = TRUE, label.size = 5, repel = TRUE) + 
  ggtitle("Cleaned T Cell UMAP (No Doublets/DN)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_umap_clean)
ggsave("Figure_T_Cell_Clean_UMAP.pdf", width = 7, height = 6)


# 2. Prepare data for proportion analysis

# Ensure Response and Timepoint columns exist (if added in sce.final previously, should be here)
# If not, let me know, I'll give you the code to add them back

# Extract metadata for plotting
meta_data <- t_cells_clean@meta.data

# Check if there are NA values (T series samples might not have Timepoint), need to exclude during analysis
meta_data_filtered <- meta_data %>% 
  filter(!is.na(Response) & !is.na(Timepoint)) %>%
  filter(Response != "Unknown") # Exclude ungrouped samples

# 3. Analysis A: Pre-treatment (Pre) - Responder vs Non-Responder

# Filter data: Pre-treatment only
data_pre <- meta_data_filtered %>% 
  filter(Timepoint == "Pre")

# Plot: 100% stacked bar chart
p_bar_pre <- ggplot(data_pre, aes(x = Response, fill = major_lineage)) +
  geom_bar(position = "fill", width = 0.7) + # position="fill" automatically calculates percentage
  scale_y_continuous(labels = scales::percent) + # Y-axis shows percentage sign
  scale_fill_manual(values = morandi_specific) + # Use Morandi colors
  theme_classic() +
  labs(title = "Baseline T Cell Composition (Pre-treatment)", 
       subtitle = "Responder vs Non-Responder",
       y = "Proportion", x = "Clinical Response") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_blank())

print(p_bar_pre)
ggsave("Figure_Proportion_Pre_R_vs_NR.pdf", width = 5, height = 4)


# 4. Analysis B: Responders - Pre vs Post treatment

# Filter data: Responders only
data_responder <- meta_data_filtered %>% 
  filter(Response == "Responder")

# Ensure X-axis order is Pre then Post
data_responder$Timepoint <- factor(data_responder$Timepoint, levels = c("Pre", "Post"))

# Plotting
p_bar_post <- ggplot(data_responder, aes(x = Timepoint, fill = major_lineage)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = morandi_specific) +
  theme_classic() +
  labs(title = "Longitudinal Changes in Responders", 
       subtitle = "Pre vs Post Treatment",
       y = "Proportion", x = "Timepoint") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_blank())

print(p_bar_post)
ggsave("Figure_Proportion_Responder_PrePost.pdf", width = 5, height = 4)

# 5. Analysis C: Non-Responders - Pre vs Post treatment

# Filter data: Non-Responders only
data_Nonresponder <- meta_data_filtered %>% 
  filter(Response == "Non-Responder")

# Ensure X-axis order is Pre then Post
data_Nonresponder$Timepoint <- factor(data_Nonresponder$Timepoint, levels = c("Pre", "Post"))

# Plotting
p_bar_npost <- ggplot(data_Nonresponder, aes(x = Timepoint, fill = major_lineage)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = morandi_specific) +
  theme_classic() +
  labs(title = "Longitudinal Changes in Non-Responders", 
       subtitle = "Pre vs Post Treatment",
       y = "Proportion", x = "Timepoint") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_blank())

print(p_bar_npost)
ggsave("Figure_Proportion_NonResponder_PrePost.pdf", width = 5, height = 4)

# ===================================================
# Calculate total T cell count before and after treatment response, related to sampling and loading, not very meaningful
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr) # Used to add significance test p-value

# Ensure using cleaned data after removal
# t_cells_clean <- ... (Object from the previous code block)


# 1. Data preparation

# Extract metadata
cell_counts_meta <- t_cells_clean@meta.data %>%
  filter(!is.na(Response) & !is.na(Timepoint) & Response != "Unknown")

# Set Morandi colors (for distinguishing groups)
group_cols <- c("Responder" = "#88A096", "Non-Responder" = "#D6B0B1") # Green vs Pink
time_cols <- c("Pre" = "#9E8E8E", "Post" = "#8F9FBC")     # Grey vs Blue

# 2. Analysis A: Pre-treatment (Pre) - Responder vs Non-Responder
# Filter pre-treatment data
data_pre_counts <- cell_counts_meta %>% filter(Timepoint == "Pre")

# --- Figure 1.1: Total T cell count (Absolute Sum) ---
p_count_total_1 <- ggplot(data_pre_counts, aes(x = Response, fill = Response)) +
  geom_bar(stat = "count", width = 0.6) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 5) + # Show specific numbers
  scale_fill_manual(values = group_cols) +
  theme_classic() +
  ylim(0, max(table(data_pre_counts$Response)) * 1.1) + # Automatically adjust Y-axis height
  labs(title = "Total T Cell Count (Pre-treatment)", 
       subtitle = "Responder vs Non-Responder (Absolute Sum)",
       y = "Number of Cells") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Figure 1.2: Average T cells per sample - Recommended! ---
# Calculate how many cells in each sample first
sample_counts_pre <- data_pre_counts %>%
  group_by(orig.ident, Response) %>%
  summarise(count = n())

p_count_boxplot_1 <- ggplot(sample_counts_pre, aes(x = Response, y = count, fill = Response)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) + # Plot points for each sample
  stat_compare_means(method = "wilcox.test", label = "p.format", vjust = -1) + # Add statistical test
  scale_fill_manual(values = group_cols) +
  theme_classic() +
  labs(title = "T Cell Infiltration per Patient (Pre)", 
       subtitle = "Does baseline infiltration predict response?",
       y = "Cells per Sample") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))


# 3. Analysis B: Responders - Pre vs Post treatment

# Filter Responder data
data_resp_counts <- cell_counts_meta %>% filter(Response == "Responder")
# Lock X-axis order
data_resp_counts$Timepoint <- factor(data_resp_counts$Timepoint, levels = c("Pre", "Post"))

# --- Figure 2.1: Total T cell count ---
p_count_total_2 <- ggplot(data_resp_counts, aes(x = Timepoint, fill = Timepoint)) +
  geom_bar(stat = "count", width = 0.6) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 5) +
  scale_fill_manual(values = time_cols) +
  theme_classic() +
  ylim(0, max(table(data_resp_counts$Timepoint)) * 1.1) +
  labs(title = "Total T Cell Expansion (Responder)", 
       subtitle = "Pre vs Post Treatment",
       y = "Number of Cells") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Figure 2.2: Average T cells per sample ---
sample_counts_resp <- data_resp_counts %>%
  group_by(orig.ident, Timepoint) %>%
  summarise(count = n())

p_count_boxplot_2 <- ggplot(sample_counts_resp, aes(x = Timepoint, y = count, fill = Timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.format", vjust = -1) + # Statistical test
  scale_fill_manual(values = time_cols) +
  theme_classic() +
  labs(title = "T Cell Infiltration Change", 
       subtitle = "Did T cells expand after treatment?",
       y = "Cells per Sample") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))


# 4. Patchwork display and save
library(patchwork)

# Patchwork: Left compares response, right compares pre/post treatment
# Top is total count, bottom is boxplot
final_plot <- (p_count_total_1 | p_count_total_2) / (p_count_boxplot_1 | p_count_boxplot_2)

print(final_plot)

ggsave("Figure_T_Cell_Counts_Statistics.pdf", width = 10, height = 8)

# ===================================================
# Found significant changes in CD4 proportion, next further split T cells into CD4 and CD8 for separate subclustering annotation
# Redefine master_palette (to prevent loss in your environment)
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
  "CD4+ Tem (Cytotoxic)"          = "#BB0021", 
  "CD4+ CTL (GNLY+)"              = "#DC0000", 
  "CD4+ Tem"                      = "#8491B4", 
  "CD4+ Proliferating"            = "#999999"
)

cols_cd8 <- c(
  # --- Naive & Stem-like ---
  "CD8+ Naive"                  = "#78C679", # Corresponds to Cluster 7
  "CD8+ Tcm (Stem-like)"        = "#41AB5D", # Corresponds to Cluster 1
  
  # --- Exhaustion-related (Cluster 8 & 0) ---
  "CD8+ Tex (GZMK+ High TOX)"   = "#D9F0A3", # Corresponds to Cluster 8 (Original Pre-Ex color)
  "CD8+ Tex (CXCL13+)"          = "#4A1486", # Corresponds to Cluster 0 (Original Terminal color)
  
  # --- Effector & Transitional (Cluster 4 & 6) ---
  "CD8+ Pre-Trm"                = "#FD8D3C", # Corresponds to Cluster 4 (Original Transitional color)
  "CD8+ Tem (GZMB+)"            = "#F16913", # Corresponds to Cluster 6
  
  # --- Resident (Cluster 2 & 5) ---
  "CD8+ Trm"                    = "#6A51A3", # Corresponds to Cluster 2 (Original Trm Exhausted color)
  "CD8+ Trm (Resting)"          = "#807DBA", # Corresponds to Cluster 5
  
  # --- Proliferating ---
  "CD8+ Proliferating"          = "#D9D9D9"  # Corresponds to Cluster 3
)

cols_other <- c(
  "Gd-T"                          = "#008B45", 
  "DN"                            = "#E1E1E1", 
  "Doublet"                       = "#E1E1E1"
)

# Combine
master_palette <- c(cols_cd4, cols_cd8, cols_other)
# Step 1: Extract and independently Re-cluster (Sub-clustering)
library(Seurat)
library(dplyr)

# Assuming your total object is called t_cells_clean (Doublets and DN removed)
# 1. Process CD8+ T cells
# Extract
cd8_obj <- subset(t_cells_clean, subset = major_lineage %in% c("CD8+", "CD8+ (Prolif)"))

# Rerun standard pipeline (Must rerun FindVariableFeatures!)
DefaultAssay(cd8_obj) <- "RNA"
cd8_obj <- NormalizeData(cd8_obj) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) # You can adjust this step, 0.4-0.8 are all fine

DimPlot(cd8_obj, reduction = "umap", label = TRUE) + ggtitle("Re-clustered CD8+T Cells")

# Annotate sub-clusters
# Define CD8 gene list
cd8_markers <- list(
  "Naive" = c("LEF1", "TCF7", "CCR7"),
  "Memory/Pre-Ex (GZMK+)" = c("GZMK", "EOMES"),
  "Progenitor Tex (Stem-like)" = c("TCF7", "PDCD1"), # TCF1+ PD1+ is the key point
  "Terminal Tex" = c("HAVCR2", "LAG3", "ENTPD1", "TOX"),
  "Tumor-Reactive (Specific)" = c("CXCL13", "LAYN"),
  "Effector (Cytotoxic)" = c("FGFBP2", "CX3CR1", "GNLY", "GZMB"),
  "Resident (Trm)" = c("ZNF683", "ITGAE"),
  "Proliferating" = c("MKI67")
)
# Plot DotPlot
DotPlot(cd8_obj, 
        features = unique(unlist(cd8_markers)), 
        group.by = "seurat_clusters", # Ensure the re-clustered IDs are used
        cols = c("#E6E6E6", "#377EB8"), # Blue tone palette
        dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "CD8+ Sub-cluster Identification", x = "Features", y = "Cluster") +
  theme(axis.text.x = element_text(face = "bold"))
# ===================================================
# Final annotation of CD8+ sub-clusters (Cluster 4 independent version)

new_cd8_ids <- c(
  "0" = "CD8+ Tex (CXCL13+)",      # Terminal exhaustion (TIM3+, CXCL13+)
  "1" = "CD8+ Tcm (Stem-like)",     # Stem-like memory (TCF7+, GZMK+)
  "2" = "CD8+ Trm",     # Exhausted resident (Hobit+, PD1+)
  "3" = "CD8+ Proliferating",       # Proliferating
  "4" = "CD8+ Pre-Trm",  # [Correction] Transitional state (GZMK+, PD1+, TCF7-)
  # Explanation: It is the bridge connecting Stem-like and Terminal, later stage than 8
  "5" = "CD8+ Trm (Resting)",       # Resting resident (Hobit+, PD1-)
  "6" = "CD8+ Tem (GZMB+)",       # Effector memory (Cytotoxic GZMB+)
  "7" = "CD8+ Naive",               # Naive
  "8" = "CD8+ Tex (GZMK+ High TOX)"  # [Relay] GZMK+, TOX high (Deepened exhaustion, loss of stemness)
)

# Rename
cd8_obj <- RenameIdents(cd8_obj, new_cd8_ids)

# Save to metadata
cd8_obj$cell_type_detailed <- Idents(cd8_obj)

# Custom style plotting
p_final <- DimPlot(cd8_obj, reduction = "umap", 
        label = TRUE, repel = TRUE, label.size = 4, 
        pt.size = 0.8, cols = cols_cd8) +
  ggtitle("CD8+ Sub-clusters (Final Annotated)")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right")

print(p_final)
ggsave("Figure_CD8_Subclusters.pdf", width = 12, height = 8)

# Or Science style (slightly brighter colors)
# DimPlot(...) + scale_color_aaas()
# Check current cell identities
print(levels(cd8_obj))


# 2. Process CD4+ T cells
# Extract
cd4_obj <- subset(t_cells_clean, subset = major_lineage == "CD4+")

# Rerun standard pipeline
DefaultAssay(cd4_obj) <- "RNA"
cd4_obj <- NormalizeData(cd4_obj) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.6) # CD4 is usually better with slightly finer clustering

DimPlot(cd4_obj, reduction = "umap", label = TRUE) + ggtitle("Re-clustered CD4+T Cells")

# Start annotation
library(Seurat)
library(ggplot2)
# 1. Prepare CD4 core Marker list
# This Panel is specifically designed to distinguish CD4 in the tumor microenvironment
cd4_panel <- list(
  "Naive / Tcm"        = c("CCR7", "LEF1", "TCF7", "SELL"),    # Naive/Central memory
  "Treg (Core)"        = c("FOXP3", "IL2RA"),                  # Regulatory T core (Must check)
  "Treg (Activated)"   = c("TNFRSF9", "TIGIT", "CTLA4"),       # Activated/Highly suppressive Treg (4-1BB, TIGIT)
  "Th1-like / Tfh"     = c("CXCL13", "BCL6", "CXCR5", "PDCD1"),# T follicular helper-like (TLS-related, star sub-cluster)
  "Th1 (Inflam)"       = c("IFNG", "TBX21", "STAT1"),          # Classic antiviral/inflammatory Th1
  "Th17 (Mucosal)"     = c("RORC", "IL17A", "KLRB1"),          # Th17 (Corresponds to previous Cluster 13)
  "CD4-CTL (Cytotoxic)"= c("GZMB", "NKG7", "GNLY"),            # CD4 with cytotoxic function
  "Proliferating"      = c("MKI67", "TOP2A")                   # Proliferating
)

# 2. Plot DotPlot
# Note: Ensure your cd4_obj has completed RunPCA -> FindClusters etc.
DotPlot(cd4_obj, 
        features = unique(unlist(cd4_panel)), 
        group.by = "seurat_clusters", # Use numerical IDs after re-clustering
        cols = c("#E6E6E6", "#D51F26"), # Grey-red palette
        dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "CD4+ T Cell Sub-cluster Identification", x = "Markers", y = "Cluster ID") +
  theme(axis.text.x = element_text(face = "bold", size = 10))


#  CD4+ Sub-cluster Annotation (Final Split Version)
new_cd4_ids_final <- c(
  # --- Treg Family ---
  "12" = "CD4+ Treg (Naive)",
  "9"  = "CD4+ Treg (Intermediate)",
  "1"  = "CD4+ Treg (Activated)",
  "11" = "CD4+ Treg (Activated)",
  
  # --- Effector/Cytotoxic Family (Separate 8 and 14) ---
  "10" = "CD4+ Tfh (CXCL13+)",
  "3"  = "CD4+ Treg (Th1-Like Transitional)",
  "6"  = "CD4+ Th17",
  
  "8"  = "CD4+ Tem (Cytotoxic)",     # [Modified] Cluster 8: Cytotoxic effector memory
  "14" = "CD4+ CTL (GNLY+)",         # [Modified] Cluster 14: Terminal CTL with strongest cytotoxicity
  
  "15" = "CD4+ Tem",                 # General effector memory
  
  # --- Naive and Interferon Family ---
  "4"  = "CD4+ Naive",
  "5"  = "CD4+ T_ISG (STAT1+)",
  "0"  = "CD4+ Tcm",
  "2"  = "CD4+ Tcm",
  "7"  = "CD4+ Tcm",
  
  # --- Proliferating ---
  "13" = "CD4+ Proliferating"
)

# Rename and save
cd4_obj <- RenameIdents(cd4_obj, new_cd4_ids_final)
cd4_obj$cell_type_detailed <- Idents(cd4_obj)

# 2. Plotting
# Plot final UMAP
p_final <- DimPlot(cd4_obj, 
                   reduction = "umap", 
                   label = TRUE, 
                   repel = TRUE, 
                   label.size = 3.5, 
                   pt.size = 0.8,    # Slightly enlarge the points for better texture
                   cols = cols_cd4) +
  ggtitle("CD4+ T Cell Sub-clusters (Final Annotated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right")

print(p_final)

# Save high-res large image
ggsave("Figure_CD4_Subclusters.pdf", width = 12, height = 8)


# 3. Merge subdivided results back into the total object (Metadata Injection)
# 3.1 Create a new column called "detailed_type", fill in major lineages first
# This way, undivided cells like Gd-T and Prolif will keep their original names and won't become NA
t_cells_clean$detailed_type <- as.character(t_cells_clean$major_lineage)

# 3.2 Use Barcode matching to fill in new CD8 names
# Logic: Find cells in cd8_obj and change their detailed_type in t_cells_clean
t_cells_clean$detailed_type[Cells(cd8_obj)] <- as.character(Idents(cd8_obj))

# 3.3 Use Barcode matching to fill in new CD4 names
t_cells_clean$detailed_type[Cells(cd4_obj)] <- as.character(Idents(cd4_obj))

# 4. Check results
# See how many detailed classifications exist now
print(table(t_cells_clean$detailed_type))

# Set default to detailed classification
Idents(t_cells_clean) <- "detailed_type"

# Plot to see the effect (Total UMAP showing sub-clusters)
p_global <- DimPlot(t_cells_clean, 
                    reduction = "umap", 
                    label = TRUE, 
                    repel = TRUE, 
                    label.size = 3,       # Global plot has many categories, make font slightly smaller to prevent overlap
                    pt.size = 0.5,        # When cell count is high, smaller points are clearer
                    cols = master_palette,      # Use the merged palette
                    raster = FALSE) +     # If cell count is extremely large (>50k), recommend setting to TRUE
  ggtitle("Global T Cell Landscape (Final Annotated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right",
        legend.text = element_text(size = 10)) # Legend font size adjustment

print(p_global)
ggsave("Figure_Total_T_cell_Annotation_UMAP.pdf", width = 15, height = 10)


# =======================================================
# 1. Organize and merge gene panel (Logical sorting version)
# High-res large DotPlot - All T cells
# We disassemble your CD4 and CD8 lists and reorder them by biological function
# This way, the resulting plot has clusters on the left and functions on the right, making it very clear

markers_clean <- unique(c(
  # --- Lineage ---
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
  
  # --- Naive & Stem-like ---
  "CCR7", "LEF1", "TCF7", "SELL", 
  
  # --- CD4 Specific (Treg, Tfh, Th17) ---
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF9",  # Treg
  "CXCL13", "BCL6", "CXCR5",                      # Tfh
  "RORC", "IL17A", "KLRB1",                       # Th17
  "TBX21", "IFNG", "STAT1",                       # Th1/ISG
  
  # --- Effector and Cytotoxic (Shared by CD4/CD8) ---
  "GZMK", "EOMES",                  # Early/Transitional effector
  "GZMB", "GNLY", "NKG7", "PRF1",   # Strong cytotoxic
  "FGFBP2", "CX3CR1",               # Terminal effector
  
  # --- Exhaustion ---
  "PDCD1", "HAVCR2", "LAG3", "ENTPD1", "TOX", "LAYN",
  
  # --- Resident ---
  "ZNF683", "ITGAE",
  
  # --- Proliferating ---
  "MKI67", "TOP2A"
))
markers_with_gdt <- unique(c(
  markers_clean,  # Previous CD4+CD8 list
  c("TRDC", "TRGC1", "TRGC2", "TRDV2") # Core markers for Gd-T
))
# Specify DotPlot order
# Extract the CD8 and CD4 name lists you defined previously
cd8_levels <- unique(new_cd8_ids)
cd4_levels <- unique(new_cd4_ids_final)
gdt_id <- "Gd-T"  # Confirm the object is indeed called this name

# Merge into the final master list (Ensure Gd-T is included!)
all_levels_ordered <- c(cd8_levels, cd4_levels, gdt_id)
# Switch identity to detailed_type first (It is a character type, contains all names)
Idents(t_cells_clean) <- "detailed_type"

# Convert it to a factor and specify the order
# This time all_levels_ordered includes "Gd-T", so it won't throw an error or generate NAs
t_cells_clean$detailed_type <- factor(t_cells_clean$detailed_type, levels = all_levels_ordered)

# Update Seurat's default identity to this ordered factor
Idents(t_cells_clean) <- t_cells_clean$detailed_type

# =======================================================
# 3. Prepare gene list (Add γδT Markers)
# =======================================================

# Ensure TRDC is in the gene list, with no duplicates
markers_final <- unique(c(
  # Previously organized CD4/CD8 global list
  markers_combined, 
  # Add Gd-T specific genes
  c("TRDC", "TRGC1", "TRGC2", "TRDV2") 
))

# =======================================================
# 4. Plotting (Perfect ending)
# =======================================================

p_final_dot <- DotPlot(t_cells_clean, 
                       features = markers_final, 
                       cols = c("#E6E6E6", "#D51F26"), 
                       dot.scale = 6) +
  RotatedAxis() +
  labs(title = "Global T Cell Landscape (Ordered & Grouped)", 
       x = "Features", y = "Cell Type") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 10), # Now the Y-axis will be strictly ordered as CD8 -> CD4 -> Gd-T
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed")
  )

print(p_final_dot)

# Provide enough height for saving to prevent text from squeezing together
ggsave("Figure_T_total_dotplot_Final_Corrected.pdf", plot = p_final_dot, width = 18, height = 12)





# =======================================================
# Save newly annotated T cell data
saveRDS(t_cells_clean, "t_cells_clean_final_annotated0310.rds")


# Plot more detailed proportion graphs
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Prepare data and color palette (Inherit from before)
meta_final <- t_cells_clean@meta.data %>%
  filter(!is.na(Response) & !is.na(Timepoint) & Response != "Unknown")
table(meta_final$detailed_type)
# Redefine master_palette (To prevent losing it in your environment)
# Fix palette (Matching names perfectly)
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
  "CD4+ Tem (Cytotoxic)"          = "#BB0021", 
  "CD4+ CTL (GNLY+)"              = "#DC0000", 
  "CD4+ Tem"                      = "#8491B4", 
  "CD4+ Proliferating"            = "#999999"
)

cols_cd8 <- c(
  # --- Naive & Stem-like ---
  "CD8+ Naive"                  = "#78C679", # Corresponds to Cluster 7
  "CD8+ Tcm (Stem-like)"        = "#41AB5D", # Corresponds to Cluster 1
  
  # --- Exhaustion-related (Cluster 8 & 0) ---
  "CD8+ Tex (GZMK+ High TOX)"   = "#D9F0A3", # Corresponds to Cluster 8 (Original Pre-Ex color)
  "CD8+ Tex (CXCL13+)"          = "#4A1486", # Corresponds to Cluster 0 (Original Terminal color)
  
  # --- Effector & Transitional (Cluster 4 & 6) ---
  "CD8+ Pre-Trm"                = "#FD8D3C", # Corresponds to Cluster 4 (Original Transitional color)
  "CD8+ Tem (GZMB+)"            = "#F16913", # Corresponds to Cluster 6
  
  # --- Resident (Cluster 2 & 5) ---
  "CD8+ Trm"                    = "#6A51A3", # Corresponds to Cluster 2 (Original Trm Exhausted color)
  "CD8+ Trm (Resting)"          = "#807DBA", # Corresponds to Cluster 5
  
  # --- Proliferating ---
  "CD8+ Proliferating"          = "#D9D9D9"  # Corresponds to Cluster 3
)
cols_other <- c(
  "Gd-T"                          = "#008B45", 
  "DN"                            = "#E1E1E1", 
  "Doublet"                       = "#E1E1E1"
)

# Combine
master_palette <- c(cols_cd4, cols_cd8, cols_other)
# 2. Self-check step (Strongly recommend running this step)
missing_colors <- setdiff(unique(meta_final$detailed_type), names(master_palette))
if(length(missing_colors) > 0) {
  cat("Warning! The following cell types have no color assigned and will appear grey or disappear in plots:\n")
  print(missing_colors)
} else {
  cat("Perfect! All cell types have been assigned colors.\n")
}

# Generate and save 5 independent plots
# ------------------------------------------------

# ==============================================================================
# Unified plotting style settings (Based on P3 mode: clean, no bold)
# ==============================================================================

# [Figure 1] Pre-treatment Global Landscape (Pre - Global)
# Comparison: Responder vs Non-Responder
data_pre <- meta_final %>% filter(Timepoint == "Pre")
p1 <- ggplot(data_pre, aes(x = Response, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 1: Baseline Landscape (All T Cells)", y = "Proportion")

print(p1)
ggsave("Fig1_Baseline_Global.pdf", plot = p1, width = 8, height = 4)


# [Figure 2] Pre-treatment CD4 Close-up (Pre - CD4 Only)
# Comparison: Responder vs Non-Responder
data_pre_cd4 <- data_pre %>% filter(major_lineage == "CD4+")
p2 <- ggplot(data_pre_cd4, aes(x = Response, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 2: Baseline CD4+ Compartment", y = "Proportion within CD4+")

print(p2)
ggsave("Fig2_Baseline_CD4_Only.pdf", plot = p2, width = 6, height = 4)


# [Figure 3] Post-treatment Global Changes (Responder - Global) --- [Baseline Style]
# Comparison: Pre vs Post
data_resp <- meta_final %>% filter(Response == "Responder")
data_resp$Timepoint <- factor(data_resp$Timepoint, levels = c("Pre", "Post")) # Lock time sequence
p3 <- ggplot(data_resp, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 3: Longitudinal Changes (Responder - All T Cells)", y = "Proportion")

print(p3)
ggsave("Fig3_Responder_Global_Change.pdf", plot = p3, width = 8, height = 4)


# [Figure 4] Post-treatment CD4 Close-up (Responder - CD4 Only)
# Comparison: Pre vs Post
data_resp_cd4 <- data_resp %>% filter(major_lineage == "CD4+")
p4 <- ggplot(data_resp_cd4, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 4: CD4+ Dynamics in Responders", y = "Proportion within CD4+")

print(p4)
ggsave("Fig4_Responder_CD4_Change.pdf", plot = p4, width = 6, height = 4)


# [Figure 5] Post-treatment CD4 Close-up (Non-Responder - CD4 Only)
# Comparison: Pre vs Post
data_nr_cd4 <- meta_final %>% 
  filter(Response == "Non-Responder") %>% 
  filter(major_lineage == "CD4+")
data_nr_cd4$Timepoint <- factor(data_nr_cd4$Timepoint, levels = c("Pre", "Post"))

p5 <- ggplot(data_nr_cd4, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 5: CD4+ Dynamics in Non-Responders", y = "Proportion within CD4+")

print(p5)
ggsave("Fig5_NonResponder_CD4_Change.pdf", plot = p5, width = 6, height = 4)


# [Figure 6] Post-treatment Global Changes (Non-Responder - Global)
# Comparison: Pre vs Post
data_nr <- meta_final %>% filter(Response == "Non-Responder")
data_nr$Timepoint <- factor(data_nr$Timepoint, levels = c("Pre", "Post"))

p6 <- ggplot(data_nr, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() + # Restore clean style, remove previous bold settings
  labs(title = "Figure 6: Longitudinal Changes (Non-Responder - All T Cells)", y = "Proportion")

print(p6)
ggsave("Fig6_NonResponder_Global_Change.pdf", plot = p6, width = 8, height = 4)


# [Figure 7] Post-treatment Global Landscape (Post - Global)
# Comparison: Responder vs Non-Responder
data_post <- meta_final %>% filter(Timepoint == "Post")
p7 <- ggplot(data_post, aes(x = Response, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() + # Restore clean style
  labs(title = "Figure 7: Post-treatment Landscape (All T Cells)", y = "Proportion")

print(p7)
ggsave("Fig7_PostTreatment_Global.pdf", plot = p7, width = 8, height = 4)

# [Figure 8] Post-treatment CD8 Close-up (Non-Responder - CD8 Only)
# Comparison: Pre vs Post
data_nr_cd8 <- meta_final %>% 
  filter(Response == "Non-Responder") %>% 
  filter(major_lineage == "CD8+")

data_nr_cd8$Timepoint <- factor(data_nr_cd8$Timepoint, levels = c("Pre", "Post"))

# Note: Object name changed to p8, title has been corrected
p8 <- ggplot(data_nr_cd8, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(title = "Figure 8: CD8+ Dynamics in Non-Responders", y = "Proportion within CD8+")

print(p8)
ggsave("Fig8_NonResponder_CD8_Change.pdf", plot = p8, width = 6, height = 4)

# [Figure 9] Post-treatment CD8 Close-up (Responder - CD8 Only)
# Purpose: Check changes in CD8 sub-clusters over time in responders
data_r_cd8 <- meta_final %>% 
  filter(Response == "Responder") %>% 
  filter(major_lineage == "CD8+")

# Ensure correct time sequence
data_r_cd8$Timepoint <- factor(data_r_cd8$Timepoint, levels = c("Pre", "Post"))

p9 <- ggplot(data_r_cd8, aes(x = Timepoint, fill = detailed_type)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = master_palette) +
  theme_classic() +
  labs(
    title = "Figure 9: CD8+ Dynamics in Responders", 
    y = "Proportion within CD8+",
    x = "Timepoint"
  )

print(p9)
ggsave("Fig9_Responder_CD8_Change.pdf", plot = p9, width = 6, height = 4)

# [Figure 10] Key Sub-cluster Aggregated Comparison: Total Trm vs Tem
# 1. Define sub-cluster list for Trm
trm_subtypes <- c("CD8+ Trm", "CD8+ Pre-Trm")

# 2. Data processing: Filter CD8 and create a new grouping column (Broad_Group)
data_grouped <- meta_final %>%
  filter(major_lineage == "CD8+") %>%
  mutate(Broad_Group = case_when(
    detailed_type %in% trm_subtypes ~ "Total Trm",  # Combine three Trms
    detailed_type == "CD8+ Tem (GZMB+)"     ~ "Tem",        # Pull out Tem separately (Ensure it's called this in your data, or is it CD8+ Teff?)
    TRUE                            ~ "Other CD8+"  # All other CD8 (Naive, Tex, etc.)
  ))

# 3. Set factor levels (Control stacking order)
# We want Trm and Tem at the bottom or the most prominent position
data_grouped$Broad_Group <- factor(data_grouped$Broad_Group, 
                                   levels = c("Total Trm", "Tem", "Other CD8+"))

data_grouped$Timepoint <- factor(data_grouped$Timepoint, levels = c("Pre", "Post"))

# 4. Define key colors (Highlight Trm with reds, Tem with blues, Other with grey)
# This makes the chart focus clear at a glance
custom_colors <- c("Total Trm" = "#E41A1C",  # Red
                   "Tem" = "#377EB8",        # Blue
                   "Other CD8+" = "grey80")  # Grey background

# 5. Plotting: Compare Responder and Non-Responder simultaneously
p10 <- ggplot(data_grouped, aes(x = Timepoint, fill = Broad_Group)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = custom_colors) +
  facet_grid(. ~ Response) + # Facet grid left/right to compare R and NR
  theme_classic() +
  labs(
    title = "Figure 10: Aggregated Dynamics: Total Trm vs Tem",
    subtitle = "Total Trm includes: Trm, and Pre-Trm",
    y = "Proportion within CD8+",
    x = "Timepoint",
    fill = "Cell Group"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "right"
  )

print(p10)
ggsave("Fig10_TotalTrm_vs_Tem_Dynamics2.pdf", plot = p10, width = 7, height = 5)

# [Figure 11] Key Sub-cluster Aggregated Comparison: CD4+ T_ISG vs Th1-Like Treg
# 1. Target sub-cluster names (Ensure exact match with strings in data)
target_isg <- "CD4+ T_ISG (STAT1+)"
target_treg <- "CD4+ Treg (Th1-Like Transitional)"

# 2. Data processing: Filter CD4 and create grouping
data_cd4_focus <- meta_final %>%
  filter(major_lineage == "CD4+") %>%
  mutate(Highlight_Group = case_when(
    detailed_type == target_isg  ~ "T_ISG (STAT1+)",      # Simplify legend name
    detailed_type == target_treg ~ "Treg (Th1-Like)",     # Simplify legend name
    TRUE                         ~ "Other CD4+"           # Other background cells
  ))

# 3. Set factor levels (Control stacking order)
# Put cells of interest at the bottom or top for easy observation
data_cd4_focus$Highlight_Group <- factor(data_cd4_focus$Highlight_Group, 
                                         levels = c("T_ISG (STAT1+)", "Treg (Th1-Like)", "Other CD4+"))

data_cd4_focus$Timepoint <- factor(data_cd4_focus$Timepoint, levels = c("Pre", "Post"))

# 4. Define colors (Use high contrast colors)
# T_ISG uses purple (representing interferon signal/inflammation), Treg uses orange/teal, background uses grey
focus_colors <- c("T_ISG (STAT1+)" = "#984EA3",   # Purple
                  "Treg (Th1-Like)" = "#FF7F00",  # Orange
                  "Other CD4+" = "grey85")        # Light grey background

# 5. Plotting
p11 <- ggplot(data_cd4_focus, aes(x = Timepoint, fill = Highlight_Group)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = focus_colors) +
  facet_grid(. ~ Response) + # Facet grid left/right comparison
  theme_classic() +
  labs(
    title = "Figure 11: CD4+ T_ISG & Th1-Like Treg Dynamics",
    subtitle = "Comparing Responders vs Non-Responders",
    y = "Proportion within CD4+",
    x = "Timepoint",
    fill = "Cell Type"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p11)
ggsave("Fig11_CD4_ISG_Treg_Dynamics.pdf", plot = p11, width = 7, height = 5)