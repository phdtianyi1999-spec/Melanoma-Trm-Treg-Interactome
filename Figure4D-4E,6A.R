#Figure4D-E
library(Seurat)
library(msigdbr)
library(dplyr)
library(pheatmap)
library(ggplot2)
getwd()
setwd("D:/ReferenceDataset/")
setwd("D:/Data_us/data/res/")
t_cells <- t_cells_clean
# ====================================================
# 1. Prepare gene sets (using msigdbr)
# ====================================================
# Get human Hallmark gene sets
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
saveRDS(h_gene_sets, file = "h_gene_sets.rds")
# Extract IFN-gamma and TGF-beta gene lists
ifng_genes <- h_gene_sets %>% 
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  pull(gene_symbol)

tgfb_genes <- h_gene_sets %>% 
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>% 
  pull(gene_symbol)

# Combine into a list. Note: the names of the list determine the general structure of the subsequent metadata column names
pathway_list <- list(
  IFNg_Score = ifng_genes,
  TGFb_Score = tgfb_genes
)

# ====================================================
# 2. Perform pathway scoring (AddModuleScore)
# ====================================================
library(Seurat)
library(dplyr)
library(pheatmap)
library(msigdbr)

# ====================================================
# 1. Ensure the complete CD4+ T cell object is used
# ====================================================
# To be safe, we extract CD4+ again from the total object t_cells
# This ensures that all cells from Responder and Non-Responder are included
all_types <- unique(t_cells$Detailed_Type)
cd4_types <- grep("^CD4\\+", all_types, value = TRUE)
cd4_all <- subset(t_cells, subset = Detailed_Type %in% cd4_types)

print(paste("Total CD4+ cells for analysis:", ncol(cd4_all)))

# ====================================================
# 2. Prepare gene sets (if not run before)
# ====================================================
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

ifng_genes <- h_gene_sets %>% 
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  pull(gene_symbol)

tgfb_genes <- h_gene_sets %>% 
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>% 
  pull(gene_symbol)

pathway_list <- list(
  IFNg_Score = ifng_genes,
  TGFb_Score = tgfb_genes
)

# ====================================================
# 3. Calculate pathway scores
# ====================================================
cd4_all <- AddModuleScore(
  object = cd4_all,
  features = pathway_list,
  name = "Global_Pathway_Score" 
)

# Rename (AddModuleScore automatically adds numeric suffixes 1 and 2)
# First clean up possible old column names to prevent confusion
existing_cols <- colnames(cd4_all@meta.data)
cols_to_remove <- c("IFNg_Response_Score", "TGFb_Signaling_Score")
cd4_all@meta.data <- cd4_all@meta.data[, !(existing_cols %in% cols_to_remove)]

# Rename newly generated columns
names(cd4_all@meta.data)[grep("Global_Pathway_Score1", names(cd4_all@meta.data))] <- "IFNg_Response_Score"
names(cd4_all@meta.data)[grep("Global_Pathway_Score2", names(cd4_all@meta.data))] <- "TGFb_Signaling_Score"

# ====================================================
# 4. Calculate means and plot
# ====================================================
# Extract data
plot_data_all <- cd4_all@meta.data %>%
  select(Detailed_Type, IFNg_Response_Score, TGFb_Signaling_Score)

# Calculate the mean for each sub-cluster
avg_scores_all <- plot_data_all %>%
  group_by(Detailed_Type) %>% 
  summarise(
    IFNg = mean(IFNg_Response_Score),
    TGFb = mean(TGFb_Signaling_Score)
  ) %>%
  as.data.frame()

# Format the data
rownames(avg_scores_all) <- avg_scores_all$Detailed_Type
avg_scores_all <- avg_scores_all[,-1]

# Plot heatmap
pheatmap(avg_scores_all,
         scale = "column",        # Scale by column: show the relative strength of each sub-cluster in the same pathway
         cluster_rows = TRUE,     # Cluster cell sub-clusters
         cluster_cols = FALSE,    # Only two columns, no clustering
         display_numbers = TRUE,  # Display numbers
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Global Pathway Activity (All CD4+ T Cells)",
         fontsize_row = 10,
         cellwidth = 45,
         cellheight = 15,
         angle_col = 45)

table(cd4_obj$Response)


# Responder Ifng Tgfb pathway activity
library(Seurat)
library(dplyr)
library(pheatmap)

# ====================================================
# 1. Filter: Keep only treatment-responsive (Responder) cells
# ====================================================
# Assuming your CD4 object is cd4_obj
# Note: Ensure the spelling of "Responder" matches exactly with your meta.data (case-sensitive)
responder_obj <- subset(cd4_all, subset = Response == "Responder")

print(paste("Number of Responder CD4 cells:", ncol(responder_obj)))

# ====================================================
# 2. Calculate the average score for each sub-cluster within the Responder group
# ====================================================
# Extract metadata
plot_data <- responder_obj@meta.data %>%
  select(Detailed_Type, IFNg_Response_Score, TGFb_Signaling_Score)

# Group by Detailed_Type and calculate the mean
avg_scores_R <- plot_data %>%
  group_by(Detailed_Type) %>% 
  summarise(
    IFNg = mean(IFNg_Response_Score),
    TGFb = mean(TGFb_Signaling_Score)
  ) %>%
  as.data.frame()

# Set row names
rownames(avg_scores_R) <- avg_scores_R$Detailed_Type
avg_scores_R <- avg_scores_R[,-1]

# ====================================================
# 3. Plot heatmap (Responder Context)
# ====================================================
pheatmap(avg_scores_R,
         scale = "column",       # Still recommend column scaling to see relative levels
         cluster_rows = T,    # Cluster to see similarities
         cluster_cols = F,   
         display_numbers = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Pathway Activity in Responders (CD4+ Subsets)",
         fontsize_row = 10,
         cellwidth = 40,
         cellheight = 15,
         angle_col = 45)

# Specify order
# 1. Define the absolute order of row names you want (force lock all Tregs together)
desired_order <- c(
  # --- Treg Family ---
  "CD4+ Treg (Activated)",
  "CD4+ Treg (Th1-Like Transitional)",
  "CD4+ Treg (Naive)",
  
  # --- Naive and Memory ---
  "CD4+ Tcm",
  "CD4+ CTL (NKG7+)",
  "CD4+ Naive",
  
  # --- Interferon response transition and terminal state ---
  "CD4+ Tcm (STAT1+)",
  "CD4+ T_ISG (STAT1+)",
  
  # --- Effector and Helper ---
  "CD4+ Tem",
  "CD4+ Tfh (CXCL13+)",
  "CD4+ Treg (Intermediate)",
  "CD4+ Th17"
)

# 2. Reorder your input matrix using this order 
# (Note: Ensure the names in desired_order match rownames(avg_scores_R) perfectly, no spelling errors)
avg_scores_R_ordered <- avg_scores_R[desired_order, ]

# 3. Rerun the heatmap code, the key is to turn off clustering
pheatmap(avg_scores_R_ordered,
         scale = "column",       
         cluster_rows = FALSE,    # [Key modification]: Turn off row clustering, strictly follow our layout
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Pathway Activity in Responders (CD4+ Subsets)",
         fontsize_row = 10,
         cellwidth = 40,
         cellheight = 15,
         angle_col = 45)

# Non-responder Ifng Tgfb pathway activity
library(Seurat)
library(dplyr)
library(pheatmap)

# ====================================================
# 1. Filter: Keep only treatment-non-responsive (Non-Responder) cells
# ====================================================
# Assuming your CD4 object is cd4_obj
non_responder_obj <- subset(cd4_all, subset = Response == "Non-Responder")

print(paste("Number of Non-Responder CD4 cells:", ncol(non_responder_obj)))

# ====================================================
# 2. Calculate the average score for each sub-cluster within the Non-Responder group
# ====================================================
plot_data_NR <- non_responder_obj@meta.data %>%
  select(Detailed_Type, IFNg_Response_Score, TGFb_Signaling_Score)

# Group by Detailed_Type and calculate the mean
avg_scores_NR <- plot_data_NR %>%
  group_by(Detailed_Type) %>% 
  summarise(
    IFNg = mean(IFNg_Response_Score),
    TGFb = mean(TGFb_Signaling_Score)
  ) %>%
  as.data.frame()

# Set row names
rownames(avg_scores_NR) <- avg_scores_NR$Detailed_Type
avg_scores_NR <- avg_scores_NR[,-1]

# ====================================================
# 3. Plot heatmap (Non-Responder Context)
# ====================================================
pheatmap(avg_scores_NR,
         scale = "column",       
         cluster_rows = TRUE,    
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Pathway Activity in Non-Responders",
         fontsize_row = 10,
         cellwidth = 40,
         cellheight = 15,
         angle_col = 45)
# Modify order
# 1. Reorder the data matrix for Non-Responders (NR)
desired_order <- c(
  # --- Treg Family ---
  "CD4+ Treg (Activated)",
  "CD4+ Treg (Th1-Like Transitional)",
  "CD4+ Treg (Naive)",
  "CD4+ T_ISG (STAT1+)",
  "CD4+ Tem",
  "CD4+ Treg (Intermediate)",
  "CD4+ Th17",
  "CD4+ CTL (NKG7+)",
  "CD4+ Naive",
  "CD4+ Tcm (STAT1+)",
  "CD4+ Tfh (CXCL13+)",
  "CD4+ Tcm"
)
avg_scores_NR_ordered <- avg_scores_NR[desired_order, ]

# 2. Plot heatmap for the NR group
pheatmap(avg_scores_NR_ordered,
         scale = "column",       
         cluster_rows = FALSE,   # [Key modification]: Turn off row clustering, strictly follow desired_order
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Pathway Activity in Non-Responders (CD4+ Subsets)", # Slightly completed the title to align with the R group
         fontsize_row = 10,
         cellwidth = 40,
         cellheight = 15,
         angle_col = 45)

#Figure6A
# 1. Load necessary packages
library(ggplot2)
library(dplyr)
library(ggrepel) # Used for more aesthetically pleasing label display, if not installed please use install.packages("ggrepel")
getwd()
setwd("D:/Data_us")
# 2. Read data
# Note: Based on your file structure, the first row is column names, the second row is metadata description, so we need to process it
data <- read.csv("drug_predictions.csv", header = TRUE)

# 3. Data cleaning
# Remove the first row (because it contains descriptive information like 'Detailed_Type...', not data)
data <- data[-1, ]

# Rename the first column to 'Drug' (it might have been empty or 'X' originally)
colnames(data)[1] <- "Drug"

# Ensure numeric columns are of numeric type (they might have been read as characters due to the presence of the second row)
data$rank <- as.numeric(as.character(data$rank))
data$logit <- as.numeric(as.character(data$logit))
data$prob <- as.numeric(as.character(data$prob))

# 4. Filter top 100 drugs
# Sort by rank, take the top 100
top100 <- data %>%
  arrange(rank) %>%
  head(100)

# Extract the top-ranked drug (Rank 0) for highlighting
top_drug <- top100[1, ]

# 5. Plotting
ggplot(top100, aes(x = prob, y = logit)) +
  # Plot scatter plot for all points
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  
  # Highlight the point for the first drug separately (red, larger)
  geom_point(data = top_drug, color = "red", size = 4) +
  
  # Label the name of the first drug
  # Use geom_label_repel to automatically avoid overlapping points, adding a box is clearer
  geom_label_repel(data = top_drug, aes(label = Drug), 
                   box.padding = 0.5, 
                   point.padding = 0.5,
                   color = "red",
                   fontface = "bold") +
  
  # Set theme and labels
  theme_minimal() +
  labs(title = "Top 100 Predicted Drugs: Activated Treg -> Th1-Like Treg",
       subtitle = "Highlighting the top-ranked candidate",
       x = "Predicted Probability (prob)",
       y = "Logit Score") +
  
  # Adjust font size
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 12))