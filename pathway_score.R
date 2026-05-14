## Change system error messages to English
Sys.setenv(LANGUAGE = "en")
## Disable conversion to factors
options(stringsAsFactors = FALSE)

library(msigdbr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(dplyr)
Sys.setenv(http_proxy = "http://127.0.0.1:8890")
Sys.setenv(https_proxy = "http://127.0.0.1:8890")
# 1. Get all human gene sets
m_df <- msigdbr(species = "Homo sapiens")
# Save as .rds file
getwd()
setwd("D:/ReferenceDataset/")
saveRDS(m_df, "msigdbr_human.rds")
setwd("D:/Data_us/")
# 2. Extract specific gene sets (build a list)
gene_lists_auto <- list()

# --- Th1 Signature (from BIOCARTA) ---
gene_lists_auto$Score_Th1_Path <- m_df %>% 
  filter(gs_name == "BIOCARTA_TH1TH2_PATHWAY") %>% 
  pull(gene_symbol)

# --- Treg Signature 
# Sticking to classic core genes is fine
gene_lists_auto$Score_Treg_Core <- c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "BATF")

# --- Hypoxia (from HALLMARK) ---
gene_lists_auto$Score_Hypoxia <- m_df %>% 
  filter(gs_name == "HALLMARK_HYPOXIA") %>% 
  pull(gene_symbol)

# --- Glycolysis (from HALLMARK) ---
gene_lists_auto$Score_Glycolysis <- m_df %>% 
  filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% 
  pull(gene_symbol)

# --- Interferon Gamma Response (from HALLMARK) ---
gene_lists_auto$Score_IFN_Response <- m_df %>% 
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  pull(gene_symbol)

# --- Apoptosis (from HALLMARK) ---
gene_lists_auto$Score_Apoptosis <- m_df %>% 
  filter(gs_name == "HALLMARK_APOPTOSIS") %>% 
  pull(gene_symbol)

# 3. Check how many genes were extracted
lapply(gene_lists_auto, length)

# 4. Run AddModuleScore
seu_sampled <- AddModuleScore(
  object = seu_sampled,
  features = gene_lists_auto,
  name = names(gene_lists_auto)
)

# --- Correct loop for fixing column names ---
# Seurat's rule: the suffix for the i-th gene set is i
for(i in seq_along(gene_lists_auto)){
  # Change the suffix here to i
  old_name <- paste0(names(gene_lists_auto)[i], i) 
  new_name <- names(gene_lists_auto)[i]
  
  # Check if the column name exists to prevent error interruption
  if(old_name %in% colnames(seu_sampled@meta.data)){
    seu_sampled[[new_name]] <- seu_sampled[[old_name]]
    seu_sampled[[old_name]] <- NULL # Delete old column
    message(paste("Renamed:", old_name, "->", new_name))
  } else {
    warning(paste("Column not found:", old_name))
  }
}

# Check results
head(seu_sampled@meta.data)

# Compare the scores of two key subclusters
target_clusters <- c("CD4+ Treg (Th1-Like Transitional)", "CD4+ T_ISG (STAT1+)")

VlnPlot(seu_sampled, 
        features = names(gene_lists_auto), 
        idents = target_clusters, 
        pt.size = 0, # Do not show points, look at the overall distribution
        ncol = 3)    # Layout

# Load necessary plotting libraries (if not already loaded)
library(ggplot2)

# Plot FeaturePlot for Th1-Like pathway score
FeaturePlot(seu_sampled,
            features = "Score_Th1_Path", # Specify the score column name to plot
            order = TRUE,                # [Key] Plot high-scoring cells on the top layer
            cols = c("#F0F0F0", "#B71C1C"), # Custom color palette: light grey to dark red to enhance contrast
            pt.size = 0.6,               # Slightly increase point size, adjust according to personal preference
            reduction = "umap"           # Default is umap, change to "tsne" if you use tsne
) +
  ggtitle("Th1-Like Pathway Score Distribution") + # Add a clear title
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Title beautification
    # axis.line = element_blank(),  # If you want a cleaner plot, you can remove axis lines and ticks
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # axis.title = element_blank()
  )

FeaturePlot(seu_sampled,
            features = "Score_Th1_Path",
            order = TRUE,
            pt.size = 0.6,
            reduction = "umap") +
  # Use magma color scale, option "A" is magma
  scale_color_viridis_c(option = "A", direction = 1) +
  ggtitle("Th1-Like Pathway Score (Magma Palette)") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(seu_sampled,
            features = "Score_Th1_Path",
            order = TRUE,         # Keep high-scoring points on the top layer
            pt.size = 0.3,
            reduction = "umap") +
  # Use ggplot2's diverging color scale function
  scale_color_gradient2(low = "#2166AC",   # Low score: dark blue
                        mid = "#F7F7F7",   # Middle: light grayish white
                        high = "#B2182B",  # High score: dark red
                        midpoint = 0) +    # Set midpoint value to 0
  ggtitle("Th1-Like Pathway Score (Diverging Scale)") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(seu_sampled,
            features = "Score_Th1_Path",
            order = TRUE,
            pt.size = 0.3,
            # Set colors: light grey -> bright red
            cols = c("lightgrey", "#E41A1C"),
            # [Key] Set cutoff values
            # "q5" and "q95" mean setting extreme values below 5% and above 95% to the minimum and maximum colors, respectively
            # This stretches the color contrast for the majority of the data in the middle
            min.cutoff = "q5",
            max.cutoff = "q95",
            reduction = "umap") +
  ggtitle("Th1-Like Pathway Score (Grey-Red with Cutoffs)") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(seu_sampled,
            features = "Score_IFN_Response",
            order = TRUE,         # Keep high-scoring points on the top layer
            pt.size = 0.3,
            reduction = "umap") +
  # Use ggplot2's diverging color scale function
  scale_color_gradient2(low = "#2166AC",   # Low score: dark blue
                        mid = "#F7F7F7",   # Middle: light grayish white
                        high = "#B2182B",  # High score: dark red
                        midpoint = 0) +    # Set midpoint value to 0
  ggtitle("Score_IFN_Response (Diverging Scale)") +
  theme(plot.title = element_text(hjust = 0.5))




# Confirm the scoring of key pathways for CD4+ Treg (Th1-Like Transitional) corresponding to IFNg reprogramming
# 1. Define gene sets
fragility_genes <- list(c("IFNG", "TBX21", "STAT1", "IL12RB1", "HIF1A", "ENTPD1"))
fragility_genes <- list(c("IFNG", "TBX21", "STAT1", "IL12RB1", "HIF1A", "ENTPD1"))
# 2. Scoring (AddModuleScore)
cd4_obj <- AddModuleScore(
  object = cd4_obj,
  features = fragility_genes,
  name = "Treg_Fragility_Score"
)
# Correct column name (Seurat automatically adds numerical suffix 1)
cd4_obj$Treg_Fragility_Score <- cd4_obj$Treg_Fragility_Score1

# 3. Extract the cell IDs for the Treg subpopulations you want to display
# Assuming your Idents are already the names in the image
target_tregs <- c("CD4+ Treg (Activated)", 
                  "CD4+ Treg (Intermediate)", 
                  "CD4+ Treg (Naive)", 
                  "CD4+ Treg (Th1-Like Transitional)")

# 4. Visualization
# Method A: Show all cells, but focus on the color intensity in the Treg area
FeaturePlot(cd4_obj, features = "Treg_Fragility_Score", 
            order = TRUE, # Let high-scoring cells float on top
            cols = c("lightgrey", "red"), 
            min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Treg Fragility Signature Activity")

# Method B (Recommended): Highlight only Treg subpopulations, other subpopulations turn grey, convenient for comparing internal heterogeneity
# split.by is not required here, but you can use it if you want to pull out Tregs separately
cells_to_plot <- WhichCells(cd4_obj, idents = target_tregs)

FeaturePlot(cd4_obj, features = "Treg_Fragility_Score", 
            cells = cells_to_plot,  # Only plot these cells
            order = TRUE,
            cols = c("lightgrey", "firebrick")) +
  NoAxes() + ggtitle("Fragility Score (Treg Subsets Only)")

# Alternative display method
library(Seurat)
library(ggplot2)
library(ggridges)
library(dplyr)

# =======================================================
# 1. Define Morandi Palette
# =======================================================
# This color palette features low saturation and grayish tones, looking elegant and soft
morandi_cols <- c(
  "CD4+ Treg (Naive)"               = "#74C476",  # Haze blue-grey
  "CD4+ Treg (Intermediate)"        = "#91D1C2",  # Warm taupe
  "CD4+ Treg (Activated)"           = "#E64B35",  # Dusty rose (corresponds to Activated, slightly redder)
  "CD4+ Treg (Th1-Like Transitional)" = "#4DBBD5"   # Mustard/Malt (corresponds to Th1-Like, slightly brighter)
)

# =======================================================
# 2. Extract and filter data
# =======================================================
# Extract scores and subpopulation information
plot_data <- FetchData(cd4_obj, vars = c("Treg_Fragility_Score", "ident")) # Be sure to check if the column name has a 1

# Modify column names for easier manipulation
colnames(plot_data) <- c("Score", "Cluster")

# [Key Step] Filter out cells with extremely low scores (cut off the sharp peak on the left)
# Based on your previous plot, the large peak on the left is around -0.5.
# We try to keep only cells with Score > -0.1 or > 0.
# This not only removes the ugly peak but also rescales the Y-axis to highlight the differences on the right.
plot_data_filtered <- plot_data %>% 
  filter(Cluster %in% names(morandi_cols)) %>% # Ensure only these 4 clusters are plotted
  filter(Score > 0.1)  # !!! This is an adjustable threshold, try 0, 0.1, or 0.2

# =======================================================
# 3. Plot the beautified ridge plot
# =======================================================
ggplot(plot_data_filtered, aes(x = Score, y = Cluster, fill = Cluster)) +
  geom_density_ridges(
    scale = 1.5,          # Let the peaks overlap slightly for a layered look
    rel_min_height = 0.01, # Cut off the excessively thin tails
    alpha = 0.8,          # Opacity, Morandi colors look better semi-transparent
    color = "white",      # Use white for peak outlines to look cleaner
    size = 0.3
  ) +
  scale_fill_manual(values = morandi_cols) + # Apply Morandi palette
  theme_classic() +
  theme(
    legend.position = "none",             # Since the Y-axis has text, the legend can be omitted
    axis.text.y = element_text(size = 12, color = "grey30"), # Y-axis font
    axis.title.x = element_text(size = 14),
    axis.line = element_line(color = "grey80"), # Lighter axis lines
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "Treg Fragility Signature (Active Cells Only)",
    x = "Fragility Score",
    y = NULL
  ) +
  # Add a dashed line representing the average to assist reading
  geom_vline(xintercept = median(plot_data_filtered$Score), 
             linetype = "dashed", color = "grey60", size = 0.5)

# 3. Plan B: Violin plot + Box plot (VlnPlot) - Most rigorous, suitable for viewing median differences
# The black line (median) in the box plot is the most telling
VlnPlot(cd4_obj, 
        features = "Treg_Fragility_Score", 
        idents = target_tregs, 
        pt.size = 0, # Do not show scatter points to avoid clutter
        cols = c("lightgrey", "lightgrey", "lightgrey", "#E41A1C")) + # Only highlight Th1-Like
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Add box plot
  ggtitle("Quantitative Comparison of Fragility")

# 4. Statistical testing 
# Extract data
score_data <- FetchData(cd4_obj, vars = c("Treg_Fragility_Score", "ident"))
# Compare only Activated and Th1-Like
sub_data <- subset(score_data, ident %in% c("CD4+ Treg (Activated)", "CD4+ Treg (Th1-Like Transitional)"))
# Wilcoxon rank-sum test
test_res <- wilcox.test(Treg_Fragility_Score ~ ident, data = sub_data)
print(paste("P-value:", test_res$p.value))