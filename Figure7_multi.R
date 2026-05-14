# ==============================================================================
# 1. Load necessary packages
# ==============================================================================
library(survival)
library(survminer)
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ==============================================================================
# 2. Data processing and merging (keep unchanged)
# ==============================================================================
# Assuming expr is a matrix and clinical is a dataframe
# 1. Ensure expr is in matrix format (required by GSVA)
expr_mat <- as.matrix(expr)

# 2. Find common samples and align
common_samples <- intersect(colnames(expr_mat), rownames(clinical))
print(paste("Number of matched samples:", length(common_samples)))

if(length(common_samples) == 0) {
  stop("Error: No samples matched! Please check if expr column names and clinical row names are consistent.")
}

expr_mat <- expr_mat[, common_samples]
clinical_matched <- clinical[common_samples, ]

expr_mat1 <- expr_mat
clinical_matched1 <- clinical_matched
rm(clinical,expr)

# A. Cohort 1 processing
expr_mat <- as.matrix(expr_mat)
colnames(expr_mat) <- paste0("C1_", colnames(expr_mat))
rownames(clinical_matched) <- paste0("C1_", rownames(clinical_matched))
clinical_matched$Batch <- "Cohort1"

# B. Cohort 2 processing
expr_mat1 <- as.matrix(expr_mat1)
colnames(expr_mat1) <- paste0("C2_", colnames(expr_mat1))
rownames(clinical_matched1) <- paste0("C2_", rownames(clinical_matched1))
clinical_matched1$Batch <- "Cohort2"

# C. Merge
common_genes <- intersect(rownames(expr_mat), rownames(expr_mat1))
expr_combined <- cbind(expr_mat[common_genes, ], expr_mat1[common_genes, ])

# D. Merge clinical information
clinical_matched <- as.data.frame(clinical_matched)
clinical_matched1 <- as.data.frame(clinical_matched1)
cols_to_keep <- c("OS", "OS.Event", "PFS", "PFS.Event", "Response", "Batch")

for(col in cols_to_keep) {
  if(!col %in% colnames(clinical_matched)) clinical_matched[[col]] <- NA
  if(!col %in% colnames(clinical_matched1)) clinical_matched1[[col]] <- NA
}

clinical_combined <- rbind(clinical_matched[, cols_to_keep], 
                           clinical_matched1[, cols_to_keep])

# ==============================================================================
# 3. Scoring (calculate three scores separately)
# ==============================================================================
trm_genes <- c("CD8A", "ITGAE", "RUNX3", "ZNF683", "GZMK", "CXCL13")
th1_treg_genes <- c("FOXP3", "IL2RA", "STAT1","HIF1A", "NR3C1",  "CD4")
combined_genes <- unique(c(trm_genes, th1_treg_genes))

calc_mean_score <- function(expr, genes) {
  valid <- intersect(rownames(expr), genes)
  if(length(valid) <= 1) return(rep(NA, ncol(expr)))
  colMeans(expr[valid, ], na.rm=TRUE)
}

# Calculate three independent scores
clinical_combined$Trm_Score <- calc_mean_score(expr_combined, trm_genes)
clinical_combined$Th1_Like_Treg_Score <- calc_mean_score(expr_combined, th1_treg_genes)
clinical_combined$Combined_Axis_Score <- calc_mean_score(expr_combined, combined_genes)

# ==============================================================================
# 4. Data cleaning and grouping
# ==============================================================================
clinical_combined$Binary_Response <- ifelse(clinical_combined$Response %in% c(1, "Responder", "CR", "PR"), 
                                            "Responder", "Non-Responder")
# ROC Outcome: 1=Responder
clinical_combined$roc_outcome <- ifelse(clinical_combined$Binary_Response == "Responder", 1, 0)

# Survival data cleaning
clinical_combined$OS <- as.numeric(as.character(clinical_combined$OS))
clinical_combined$OS.Event <- as.numeric(as.character(clinical_combined$OS.Event))
clinical_combined$PFS <- as.numeric(as.character(clinical_combined$PFS))
clinical_combined$PFS.Event <- as.numeric(as.character(clinical_combined$PFS.Event))

# Median grouping (for survival analysis)
# Here use the median of "Combined_Axis_Score"
med_val <- median(clinical_combined$Combined_Axis_Score, na.rm = TRUE)
clinical_combined$Combined_Group <- ifelse(clinical_combined$Combined_Axis_Score > med_val, "High", "Low")
clinical_combined$Combined_Group <- factor(clinical_combined$Combined_Group, levels = c("High", "Low"))

# ==============================================================================
# 5. [Core modification] Define analysis function (separate cohorts: violin plot + multiple ROC)
# ==============================================================================

analyze_cohort_full <- function(data, cohort_name) {
  
  # Extract subset with response data
  plot_data <- subset(data, !is.na(Binary_Response))
  
  # --- Part A: Violin Plot - Three-in-one display ---
  # Define internal plotting function
  my_violin <- function(y_col, title_sub) {
    ggviolin(plot_data, x = "Binary_Response", y = y_col,
             fill = "Binary_Response", palette = c("#E7B800", "#2E9FDF"),
             add = "boxplot", add.params = list(fill = "white", width = 0.1), # Nested boxplot
             title = title_sub, xlab = "Response", ylab = "Score") +
      stat_compare_means(method = "wilcox.test", label.x.npc = "center", label.y.npc = "bottom")
  }
  
  # Plot three graphs separately
  p1 <- my_violin("Trm_Score", "Trm Score")
  p2 <- my_violin("Th1_Like_Treg_Score", "Th1-Treg Score")
  p3 <- my_violin("Combined_Axis_Score", "Combined Score")
  
  # Patchwork display (one row, three graphs) - This is a very beautiful Figure panel
  print(ggarrange(p1, p2, p3, ncol = 3, nrow = 1, 
                  common.legend = TRUE, legend = "top",
                  labels = c("A", "B", "C"))) # Automatically add A, B, C labels
  
  
  # --- Part B: Multiple ROC analysis (Comparison) ---
  if(length(unique(plot_data$roc_outcome)) == 2) {
    
    # Calculate three ROC objects
    roc_trm <- roc(plot_data$roc_outcome, plot_data$Trm_Score, quiet=TRUE)
    roc_treg <- roc(plot_data$roc_outcome, plot_data$Th1_Like_Treg_Score, quiet=TRUE)
    roc_comb <- roc(plot_data$roc_outcome, plot_data$Combined_Axis_Score, quiet=TRUE)
    
    # Extract AUC values
    auc_trm <- round(auc(roc_trm), 3)
    auc_treg <- round(auc(roc_treg), 3)
    auc_comb <- round(auc(roc_comb), 3)
    
    # Plotting
    # 1. Draw the first line (Combined - thick purple line)
    plot(roc_comb, col="purple", lwd=3, legacy.axes=TRUE, grid=TRUE,
         main = paste0(cohort_name, ": ROC Analysis Comparison"))
    
    # 2. Add the second line (Trm - red dashed line)
    plot(roc_trm, col="red", lwd=2, add=TRUE, lty=2)
    
    # 3. Add the third line (Th1-Treg - blue dashed line)
    plot(roc_treg, col="blue", lwd=2, add=TRUE, lty=2)
    
    # 4. Add legend (display AUC)
    legend("bottomright", 
           legend = c(paste0("Combined (AUC=", auc_comb, ")"),
                      paste0("Trm (AUC=", auc_trm, ")"),
                      paste0("Th1-Treg (AUC=", auc_treg, ")")),
           col = c("purple", "red", "blue"),
           lwd = c(3, 2, 2),
           lty = c(1, 2, 2),
           bty = "n") # Remove legend border
    
  } else {
    print(paste("Skipping ROC for", cohort_name, "- Not enough classes."))
  }
}

# ==============================================================================
# 6. Execute analysis (separate cohorts)
# ==============================================================================
data_c1 <- subset(clinical_combined, Batch == "Cohort1")
data_c2 <- subset(clinical_combined, Batch == "Cohort2")

print("========== Analyzing Cohort 1 ==========")
analyze_cohort_full(data_c1, "Cohort 1")

print("========== Analyzing Cohort 2 ==========")
analyze_cohort_full(data_c2, "Cohort 2")


# ==============================================================================
# 7. Survival analysis (keep picture mode - combined cohorts have the strongest statistical power)
# ==============================================================================
plot_fancy_survival <- function(data, time_col, event_col, title_text) {
  
  # Build formula
  f <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ Combined_Group"))
  
  # Fitting
  fit <- survfit(f, data = data)
  
  # Inject formula (key to fixing ggsurvplot error)
  fit$call$formula <- f
  fit$call$data <- data
  
  # Plotting (highly customized)
  ggsurvplot(fit, data = data,
             pval = TRUE,             # Show P value
             pval.method = TRUE,      # Show Log-rank text
             conf.int = F,         # Show confidence interval shadow
             risk.table = TRUE,       # Show risk table
             risk.table.col = "strata", # Risk table color
             risk.table.height = 0.25,
             palette = c("#E7B800", "#2E9FDF"), # Classic color palette
             surv.median.line = "hv", # Median survival line
             ggtheme = theme_classic(),
             title = title_text,
             xlab = "Time (Months)",
             legend.labs = c("High Score", "Low Score"))
}

# Plot survival analysis of combined cohorts
print(plot_fancy_survival(clinical_combined, "OS", "OS.Event", "Overall Survival (Combined Cohorts)"))
print(plot_fancy_survival(clinical_combined, "PFS", "PFS.Event", "Progression-Free Survival (Combined Cohorts)"))


# Combine signature scoring
# ==============================================================================
# 1. Load necessary packages
# ==============================================================================
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ==============================================================================
# 2. Data preparation (if you have already run data merging before, this part will be fast)
# ==============================================================================
# To ensure the code runs independently, necessary data processing steps are included here
# ------------------------------------------------------------------------------
# A. Cohort 1 processing
expr_mat_c1 <- as.matrix(expr_mat)
colnames(expr_mat_c1) <- paste0("C1_", colnames(expr_mat)) # Prevent ID conflicts
rownames(clinical_matched) <- paste0("C1_", rownames(clinical_matched))
clinical_matched$Batch <- "Cohort1"

# B. Cohort 2 processing
expr_mat_c2 <- as.matrix(expr_mat1)
colnames(expr_mat_c2) <- paste0("C2_", colnames(expr_mat1))
rownames(clinical_matched1) <- paste0("C2_", rownames(clinical_matched1))
clinical_matched1$Batch <- "Cohort2"

# C. Merge
common_genes <- intersect(rownames(expr_mat_c1), rownames(expr_mat_c2))
expr_combined <- cbind(expr_mat_c1[common_genes, ], expr_mat_c2[common_genes, ])

# D. Merge clinical information
clinical_matched <- as.data.frame(clinical_matched)
clinical_matched1 <- as.data.frame(clinical_matched1)
cols_to_keep <- c("Response", "Batch") # Only need these two columns

# Complete columns
for(col in cols_to_keep) {
  if(!col %in% colnames(clinical_matched)) clinical_matched[[col]] <- NA
  if(!col %in% colnames(clinical_matched1)) clinical_matched1[[col]] <- NA
}

clinical_combined <- rbind(clinical_matched[, cols_to_keep], 
                           clinical_matched1[, cols_to_keep])

# ==============================================================================
# 3. Only calculate Combined Signature score (average expression)
# ==============================================================================
trm_genes <- c("CD8A", "ITGAE", "CD69", "RUNX3", "ZNF683", "PRDM1", "CXCL13")
th1_treg_genes <- c("FOXP3", "IL2RA", "TBX21", "CXCR3", "CTLA4", "LAG3")
combined_genes <- unique(c(trm_genes, th1_treg_genes))

# Scoring function
calc_mean_score <- function(expr, genes) {
  valid <- intersect(rownames(expr), genes)
  if(length(valid) <= 1) return(rep(NA, ncol(expr)))
  colMeans(expr[valid, ], na.rm=TRUE)
}

# Calculate Combined Score
clinical_combined$Combined_Axis_Score <- calc_mean_score(expr_combined, combined_genes)

# ==============================================================================
# 4. Data cleaning (standardize Response)
# ==============================================================================
# Define Responder
clinical_combined$Binary_Response <- ifelse(clinical_combined$Response %in% c(1, "Responder", "CR", "PR"), 
                                            "Responder", "Non-Responder")
# ROC needs 0/1 (1=Responder)
clinical_combined$roc_outcome <- ifelse(clinical_combined$Binary_Response == "Responder", 1, 0)

# ==============================================================================
# 5. [Core] Plot cohorts separately (only plot Combined Score)
# ==============================================================================

# Define dedicated plotting function
plot_combined_only <- function(data, cohort_name) {
  
  # Extract valid data
  plot_data <- subset(data, !is.na(Binary_Response))
  
  # --- A. Violin Plot ---
  p_violin <- ggviolin(plot_data, x = "Binary_Response", y = "Combined_Axis_Score",
                       fill = "Binary_Response", 
                       palette = c("#E7B800", "#2E9FDF"), # Classic yellow-blue palette
                       add = "boxplot", 
                       add.params = list(fill = "white", width = 0.1), # Nested white boxplot
                       title = paste0(cohort_name, ": Combined Signature Score"),
                       xlab = "Response Group", ylab = "Mean Expression Score") +
    stat_compare_means(method = "wilcox.test", 
                       label.x.npc = "center", 
                       label.y.npc = "bottom") # P-value displayed at the bottom
  
  print(p_violin)
  
  # --- B. ROC Curve ---
  if(length(unique(plot_data$roc_outcome)) == 2) {
    # Calculate ROC
    roc_obj <- roc(plot_data$roc_outcome, plot_data$Combined_Axis_Score, quiet=TRUE)
    auc_val <- round(auc(roc_obj), 3)
    
    # Plotting
    plot(roc_obj, 
         col = "purple", lwd = 3,       # Thick purple line
         legacy.axes = TRUE, grid = TRUE,
         main = paste0(cohort_name, ": ROC Curve (Combined Signature)"))
    
    # Add legend
    legend("bottomright", 
           legend = paste0("AUC = ", auc_val),
           col = "purple", lwd = 3, 
           bty = "n") # No border
    
  } else {
    print(paste("Skipping ROC for", cohort_name, "- Not enough classes (need both R and NR)."))
  }
}

# ==============================================================================
# 6. Execute plotting
# ==============================================================================
data_c1 <- subset(clinical_combined, Batch == "Cohort1")
data_c2 <- subset(clinical_combined, Batch == "Cohort2")

print("========== Plotting Cohort 1 ==========")
plot_combined_only(data_c1, "Cohort 1")

print("========== Plotting Cohort 2 ==========")
plot_combined_only(data_c2, "Cohort 2")