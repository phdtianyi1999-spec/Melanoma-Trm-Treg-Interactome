# ==============================================================================
# 1. Load necessary packages
# ==============================================================================
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ==============================================================================
# 2. Data cleaning and alignment (Key modification: remove NA)
# ==============================================================================
# Ensure correct format
exp_mat <- as.matrix(exp)
clinical <- as.data.frame(sample)

# [New Step] Remove samples where Response is NA
print(paste("Number of samples before cleaning:", nrow(clinical)))
clinical <- clinical[!is.na(clinical$Response), ]
print(paste("Number of samples after removing Response=NA:", nrow(clinical)))

# Extract sample names (assuming Sample column is the sample ID)
rownames(clinical) <- clinical$Sample

# Find common samples (take intersection)
common_samples <- intersect(colnames(exp_mat), rownames(clinical))
print(paste("Final number of matched common samples:", length(common_samples)))

if(length(common_samples) == 0) {
  stop("Error: No samples matched! Please check if exp column names and sample$Sample are consistent.")
}

# Align data
exp_sub <- exp_mat[, common_samples]
clin_sub <- clinical[common_samples, ]

# ==============================================================================
# 3. Calculate Combined Signature Score (Average Expression)
# ==============================================================================
# Define gene sets
trm_genes <- c("CD8A", "ITGAE", "CD69", "RUNX3", "ZNF683", "PRDM1", "CXCL13")
th1_treg_genes <- c("FOXP3", "IL2RA", "TBX21", "CXCR3", "CTLA4", "LAG3")
combined_genes <- unique(c(trm_genes, th1_treg_genes))

# Scoring function (kept unchanged)
calc_mean_score <- function(expression_matrix, gene_list) {
  valid_genes <- intersect(rownames(expression_matrix), gene_list)
  if(length(valid_genes) <= 1) return(rep(NA, ncol(expression_matrix)))
  return(colMeans(expression_matrix[valid_genes, ], na.rm = TRUE))
}

# 1. Calculate original score
clin_sub$Combined_Axis_Score <- calc_mean_score(exp_sub, combined_genes)

# --------------------------------------------------------------------------
# [New] Scale scores to between -2 and 2
# --------------------------------------------------------------------------
# Define normalization function
scale_to_range <- function(x, new_min = -2, new_max = 2) {
  # Prevent division by 0 (if all values are the same)
  if(max(x, na.rm=TRUE) == min(x, na.rm=TRUE)) { return(rep(0, length(x))) }
  
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) * (new_max - new_min) + new_min
}

# Generate new scaled column (recommended to keep original column and create a new one for plotting)
clin_sub$Scaled_Score <- scale_to_range(clin_sub$Combined_Axis_Score, -2, 2)

# Check the results
print("Score range after scaling:")
print(range(clin_sub$Scaled_Score, na.rm = TRUE))
# ==============================================================================
# 4. Define groups (PRCR vs SD+PD)
# ==============================================================================
# Based on your description: Response contains PD, PRCR, SD
# Responder = PRCR
# Non-Responder = PD, SD

clin_sub$Binary_Response <- ifelse(clin_sub$Response %in% c("PRCR", "CR", "PR"), 
                                   "Responder", 
                                   "Non-Responder")

# Check grouping status again after removing NAs
print("Final grouping statistics:")
print(table(clin_sub$Binary_Response))

# Define ROC variable (1=Responder)
clin_sub$roc_outcome <- ifelse(clin_sub$Binary_Response == "Responder", 1, 0)

# ==============================================================================
# 5. Plotting
# ==============================================================================

# --- A. Violin Plot ---
# Ensure there are scores and groupings
plot_data <- subset(clin_sub, !is.na(Scaled_Score)) # Note: use Scaled_Score to check for NA here

p_violin <- ggviolin(plot_data, 
                     x = "Binary_Response", 
                     y = "Scaled_Score",  # [Modified] Change this to the scaled column name
                     fill = "Binary_Response", 
                     palette = c("#2E9FDF","#E7B800"), 
                     add = "boxplot", 
                     add.params = list(fill = "white", width = 0.1), 
                     title = "Combined Signature Score (Scaled -2 to 2)", # Change the title
                     xlab = "Response Group", 
                     ylab = "Scaled Expression Score (-2 to 2)") + # Change the Y-axis label
  stat_compare_means(method = "wilcox.test", 
                     label.x.npc = "center", 
                     label.y.npc = "bottom")

print(p_violin)

# --- B. ROC Curve (Recommended to still use original scores, or scaled scores are fine, AUC value is the same) ---
if(length(unique(plot_data$roc_outcome)) == 2) {
  # Note: Linear transformation does not change the AUC value and shape of ROC, either is fine
  roc_obj <- roc(plot_data$roc_outcome, plot_data$Scaled_Score, quiet=TRUE) 
  auc_val <- round(auc(roc_obj), 3)
  
  plot(roc_obj, 
       col = "purple", lwd = 3, legacy.axes = TRUE, grid = TRUE,
       main = "ROC Curve (New Cohort)")
  
  legend("bottomright", 
         legend = paste0("AUC = ", auc_val),
         col = "purple", lwd = 3, 
         bty = "n")
  
} else {
  print("Notice: Insufficient sample classification (must contain both Responder and Non-Responder) to plot ROC.")
}