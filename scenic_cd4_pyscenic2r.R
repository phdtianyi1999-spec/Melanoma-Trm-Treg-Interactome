library(SCopeLoomR)
library(SCENIC)
library(Seurat)
getwd()
cd4_obj <- seu_sampled
# 1. Load the loom file
loom <- open_loom("cd4_scenic_output.loom")

# 2. Extract Regulon activity score matrix (AUC)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
library(AUCell)  # Ensure this package is loaded

# 1. Extract the actual numeric matrix from the aucellResults object
# getAUC is a dedicated function provided by the AUCell package
auc_matrix <- getAUC(regulonAUC)

# 2. Check the dimensions of the extracted matrix
# Should be 105 (Regulons) x 8896 (Cells)
print(dim(auc_matrix))

# 3. Clean row names (remove suffix)
rownames(auc_matrix) <- gsub(" \\(.*", "", rownames(auc_matrix))

# 4. Store directly into Seurat (do not transpose!)
# Seurat requires: rows are features (Regulons), columns are cells. Your data is already in this format.
cd4_obj[["SCENIC"]] <- CreateAssayObject(data = auc_matrix)

# 5. Verify if successful
print(cd4_obj[["SCENIC"]])

# 6. Visualization
DefaultAssay(cd4_obj) <- "SCENIC"
FeaturePlot(cd4_obj, features = c("BATF", "STAT1", "FOXP3"), reduction = "umap")

# Plot heatmap: View the most highly activated TFs within each CD4 subpopulation
# 1. Switch to the SCENIC assay (assuming you have already stored it)
DefaultAssay(cd4_obj) <- "SCENIC"
all_regulons <- rownames(cd4_obj)
cd4_obj <- ScaleData(cd4_obj, features = all_regulons)
# 2. Find Marker Regulons for each subpopulation
# Note: Set logfc.threshold a bit smaller because the magnitude of change in AUC scores is not as large as gene expression
scenic_markers <- FindAllMarkers(cd4_obj, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.05)

# 3. Select the top 5-10 most significant TFs for each subpopulation
top_regulons <- scenic_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

# 1. Define your custom color palette
cols_cd4 <- c(
  "CD4+ Naive"                    = "#74C476", 
  "CD4+ Tcm"                      = "#3C5488", 
  "CD4+ T_ISG (STAT1+)"           = "#631879", 
  "CD4+ Treg (Naive)"             = "#00A087", 
  "CD4+ Treg (Intermediate)"      = "#91D1C2", 
  "CD4+ Treg (Activated)"         = "#E64B35", 
  "CD4+ Tfh (CXCL13+)"            = "#B09C85", 
  "CD4+ Th17"                     = "#F39B7F", 
  "CD4+ Treg (Th1-Like Transitional)" = "#4DBBD5", 
  "CD4+ Tem (Cytotoxic)"          = "#BB0021", 
  "CD4+ CTL (GNLY+)"              = "#DC0000", 
  "CD4+ Tem"                      = "#8491B4", 
  "CD4+ Proliferating"            = "#999999"
)

# 2. Plot the heatmap
DoHeatmap(cd4_obj, 
          features = top_regulons$gene, # Ensure these are the TF names you want to display
          group.by = "ident",
          group.colors = cols_cd4,      # Your custom color palette
          label = FALSE,                # Remove top text
          size = 0,                     # Remove gaps
          raster = FALSE) +             # Vector graphic for higher clarity
  
  # --- Key: Black-Gold color scheme ---
  # Low values are black, median values are gold, high values are bright yellow
  scale_fill_gradientn(colors = c("black", "#B8860B", "#FFD700", "#FFFFE0"),
                       name = "Activity") +
  
  theme(
    # Y-axis font (TF names)
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    
    # Legend font
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    
    # Top title
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    
    # Remove X-axis clutter
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(title = "Top Regulon Activity (Black-Gold Heatmap)")

library(SCENIC)
library(AUCell)
library(ggplot2)
library(ggrepel)

# 1. Extract AUC matrix (Rows=Regulons, Columns=Cells)
auc_mtx <- GetAssayData(cd4_obj, assay = "SCENIC", slot = "data")

# 2. Get cell classification information (Idents)
cell_clusters <- setNames(as.character(Idents(cd4_obj)), names(Idents(cd4_obj)))

# 3. Calculate RSS (Regulon Specificity Score)
# This step calculates the specificity score of each TF in each subpopulation
rss <- calcRSS(AUC = auc_mtx, cellAnnotation = cell_clusters)

# 4. Check the rss object, it is a matrix (Rows=Regulons, Columns=Clusters)
head(rss)
# Define a plotting function to perfectly replicate the style of your uploaded image
# --- Optimized plotting function ---
plot_rss_rank_distributed <- function(rss_mat, cluster_name, top_n = 5) {
  
  # 1. Organize data
  df <- data.frame(
    Regulon = rownames(rss_mat),
    Score = rss_mat[, cluster_name]
  )
  df <- df[order(df$Score, decreasing = TRUE), ]
  df$Rank <- 1:nrow(df)
  
  # 2. Mark Top N
  df$Label <- NA
  df$Label[1:top_n] <- df$Regulon[1:top_n]
  
  # Calculate maximum value for X-axis
  max_rank <- max(df$Rank)
  
  # 3. Plotting
  ggplot(df, aes(x = Rank, y = Score)) +
    # Background points: Blue, moderate size
    geom_point(color = "#0072B5", alpha = 0.5, size = 3) + 
    
    # Highlighted points: Red, slightly larger
    geom_point(data = df[1:top_n, ], aes(x = Rank, y = Score), 
               color = "#BC3C29", size = 5) +
    
    # --- Label settings (Core modification) ---
    geom_text_repel(aes(label = Label),
                    data = df[1:top_n, ],
                    
                    # [Key 1] Remove nudge_x and direction parameters
                    # Let ggrepel freely find directions (up, down, left, right)
                    
                    # [Key 2] Increase padding distance to encourage labels to spread out
                    box.padding = 1.0,    # Padding around the text box; larger means more spread out
                    point.padding = 0.6,  # Distance between text and data points
                    
                    # Ensure all labels are displayed
                    max.overlaps = Inf,
                    
                    # Line settings: Always show segments
                    min.segment.length = 0, 
                    segment.size = 0.5,
                    segment.color = "grey60",
                    
                    size = 5, # Font size
                    fontface = "bold",
                    color = "black") +
    
    # Slightly reduce X-axis range since labels are no longer forced to the right
    scale_x_continuous(limits = c(0, max_rank * 1.1)) +
    
    labs(title = paste0(cluster_name, "\nTop ", top_n, " Regulons"),
         y = "Specificity Score (RSS)", x = "Regulons Rank") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      # Increase plot margins slightly to prevent spread labels from being cut off
      plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")
    )
}
# Display plots
plot_list[[1]] # View first plot
plot_list[[2]] # View second plot
plot_list[[3]] # View third plot

# Or arrange them together
library(patchwork)
wrap_plots(plot_list, ncol = 3)



# Plot scatter plots of AUC activity rankings
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- [1] Calculate average AUC activity for each subpopulation ---
# --- [1] Corrected calculation function (prevents names from being changed) ---
calc_mean_auc <- function(seurat_obj, assay = "SCENIC") {
  # Extract AUC matrix
  auc_mtx <- GetAssayData(seurat_obj, assay = assay, layer = "data")
  if(nrow(auc_mtx) == 0) auc_mtx <- GetAssayData(seurat_obj, assay = assay, slot = "data")
  
  # Get cell grouping
  idents <- Idents(seurat_obj)
  
  # Calculate the average value for each Cluster
  # [Key modification] Add check.names = FALSE in data.frame
  mean_auc_df <- data.frame(t(apply(auc_mtx, 1, function(x) tapply(x, idents, mean))), 
                            check.names = FALSE)
  
  return(mean_auc_df)
}

# --- [2] Recalculate matrix ---
activity_matrix <- calc_mean_auc(cd4_obj)

# Check if the column names are correct (should contain spaces and plus signs now)
print(colnames(activity_matrix))

# Execute calculation
activity_matrix <- calc_mean_auc(cd4_obj)


# --- [2] Plotting function (reuse previous beautiful layout) ---
plot_activity_rank <- function(activity_mat, cluster_name, top_n = 5) {
  
  # Organize data
  # Note: Here we are taking Mean AUC, not RSS
  df <- data.frame(
    Regulon = rownames(activity_mat),
    Score = activity_mat[[cluster_name]] # Get the corresponding column
  )
  
  # Sort (from high to low)
  df <- df[order(df$Score, decreasing = TRUE), ]
  df$Rank <- 1:nrow(df)
  
  # Mark Top N
  df$Label <- NA
  df$Label[1:top_n] <- df$Regulon[1:top_n]
  
  # Plotting
  ggplot(df, aes(x = Rank, y = Score)) +
    # Background points: Use purple tones to differentiate
    geom_point(color = "#7E6148", alpha = 0.5, size = 3) + 
    
    # Highlighted points: Dark purple
    geom_point(data = df[1:top_n, ], aes(x = Rank, y = Score), 
               color = "#E64B35", size = 5) +
    
    # Labels (smart repel)
    geom_text_repel(aes(label = Label),
                    data = df[1:top_n, ],
                    box.padding = 1.0, 
                    point.padding = 0.5,
                    min.segment.length = 0,
                    segment.color = "grey50",
                    max.overlaps = Inf,
                    size = 5, fontface = "bold", color = "black") +
    
    # Axis settings
    scale_x_continuous(limits = c(0, max(df$Rank) * 1.1)) +
    
    # Title
    labs(title = paste0(cluster_name, "\nTop Active Regulons (Mean AUC)"),
         y = "Mean AUC Activity", x = "Regulon Rank") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}

# --- [3] Plot the three clusters you want to see ---

# 1. Th1-Like Treg
p_act_1 <- plot_activity_rank(activity_matrix, "CD4+ Treg (Th1-Like Transitional)", top_n = 5)

# 2. T_ISG
p_act_2 <- plot_activity_rank(activity_matrix, "CD4+ T_ISG (STAT1+)", top_n = 5)

# 3. Activated Treg
p_act_3 <- plot_activity_rank(activity_matrix, "CD4+ Treg (Activated)", top_n = 5)

# Display image
print(p_act_1)

# Recommended operation: Plot an "RSS vs AUC" scatter plot
library(ggplot2)
library(ggrepel)

# 1. Prepare data
# Assume target_cluster is "CD4+ Treg (Th1-Like Transitional)"
# rss is the RSS matrix just calculated
# activity_matrix is the Mean AUC matrix just calculated

df_scatter <- data.frame(
  Regulon = rownames(rss),
  RSS = rss[, "CD4+ Treg (Th1-Like Transitional)"],
  AUC = activity_matrix[, "CD4+ Treg (Th1-Like Transitional)"]
)

# 2. Mark the genes you want to highlight (Union of Top RSS and Top AUC)
top_rss_genes <- rownames(df_scatter)[order(df_scatter$RSS, decreasing = T)][1:5]
top_auc_genes <- rownames(df_scatter)[order(df_scatter$AUC, decreasing = T)][1:5]
# Plus your biologically relevant genes of interest
manual_genes <- c("TBX21(+)", "STAT1(+)", "HIF1A(+)", "FOXP3(+)")

label_genes <- unique(c(top_rss_genes, top_auc_genes, manual_genes))
df_scatter$Label <- ifelse(df_scatter$Regulon %in% label_genes, df_scatter$Regulon, NA)

# 3. Plotting
ggplot(df_scatter, aes(x = AUC, y = RSS)) +
  geom_point(color = "grey80", alpha = 0.6) +
  
  # Highlight Top genes
  geom_point(data = subset(df_scatter, !is.na(Label)), 
             color = "red", size = 3) +
  
  geom_text_repel(aes(label = Label), max.overlaps = Inf) +
  
  labs(title = "Specificity (RSS) vs Activity (AUC)",
       x = "Mean Activity (AUC)", 
       y = "Specificity (RSS)") +
  theme_classic()

# View the highest expressed regulons in the two subpopulations of interest
# ==============================================================================
# 0. Load necessary packages
# ==============================================================================
library(Seurat)
library(AUCell) 
library(dplyr)

# ==============================================================================
# 1. Set targets
# ==============================================================================
seurat_obj <- cd4_obj 

# Your target subpopulations (names must exactly match those in Idents(seurat_obj))
target_clusters <- c("CD4+ Treg (Th1-Like Transitional)", "CD4+ T_ISG (STAT1+)","CD4+ Treg (Naive)","CD4+ Treg (Activated)")
top_n_cutoff <- 50  # Extract the top 50

# ==============================================================================
# 2. Prepare data: Calculate mean activity (Mean AUC)
# ==============================================================================
message("Extracting AUC matrix...")
# Extract AUC matrix (ensure the assay name is correct, usually "SCENIC" or "AUC")
auc_mtx <- GetAssayData(seurat_obj, assay = "SCENIC", layer = "data")

message("Calculating Mean AUC...")
# Calculate the average AUC for each subpopulation
# check.names = FALSE is crucial to prevent R from changing "CD4+ Treg..." to "CD4..Treg..."
mean_auc_df <- data.frame(t(apply(auc_mtx, 1, function(x) {
  tapply(x, Idents(seurat_obj), mean)
})), check.names = FALSE)

# ==============================================================================
# 3. [Key fix] Calculate RSS matrix on the fly
# ==============================================================================
message("Calculating RSS (Regulon Specificity Score)... this may take a few seconds...")

# Prepare various inputs needed for RSS
rss_auc <- auc_mtx
rss_clusters <- Idents(seurat_obj)

# Run calcRSS
rss_result <- calcRSS(AUC = rss_auc, cellAnnotation = rss_clusters)

# Ensure it is in standard matrix format
rss_matrix <- as.matrix(rss_result)

# ***Key check***: Print the column names of the RSS matrix to confirm your target subpopulations are included
message("RSS calculation complete. Check if column names contain target subpopulations:")
print(intersect(target_clusters, colnames(rss_matrix)))

# ==============================================================================
# 4. Core loop: Extract Top N RSS & Top N Activity and take the union
# ==============================================================================

# Create a list to store results, convenient for plotting heatmaps later
union_regulons_list <- list()

cat("\n")
cat("##################################################################\n")
cat("      🚀 Key Regulon Extraction Report (Top", top_n_cutoff, "Union)\n")
cat("##################################################################\n\n")

for (cluster in target_clusters) {
  
  # --- A. Find Top Activity (based on Mean AUC) ---
  # Find the column in mean_auc_df, sort in descending order, take top N
  if (cluster %in% colnames(mean_auc_df)) {
    top_activity <- rownames(mean_auc_df)[order(mean_auc_df[[cluster]], decreasing = TRUE)][1:top_n_cutoff]
  } else {
    warning(paste("Subpopulation not found in Mean AUC:", cluster))
    top_activity <- c()
  }
  
  # --- B. Find Top Specificity (based on RSS scores) ---
  # Now rss_matrix is freshly calculated, the names will definitely match
  if (cluster %in% colnames(rss_matrix)) {
    top_rss <- rownames(rss_matrix)[order(rss_matrix[, cluster], decreasing = TRUE)][1:top_n_cutoff]
  } else {
    warning(paste("Subpopulation not found in RSS matrix (please check punctuation):", cluster))
    top_rss <- c()
  }
  
  # --- C. Take the union ---
  union_set <- unique(c(top_activity, top_rss))
  union_regulons_list[[cluster]] <- union_set
  
  # --- D. Print results ---
  cat(" Subpopulation:", cluster, "\n")
  cat("------------------------------------------------------------------\n")
  
  cat(" [Top Activity] (High activity): \n")
  cat(paste(top_activity, collapse = ", "), "\n\n")
  
  cat(" [Top RSS] (High specificity - RSS): \n")
  cat(paste(top_rss, collapse = ", "), "\n\n")
  
  cat("[Final Union] (Total", length(union_set), "): \n")
  cat(paste(union_set, collapse = ", "), "\n")
  cat("\n\n")
}
cat("##################################################################\n")
cat("✅ Extraction complete!\n")



getwd()
#===================================================
# Generate table from results
# ================= STEP 1: Prepare containers =================
# Use a List to temporarily store the results for each Cluster, as Lists allow different lengths
all_regulons_list <- list()

# ================= STEP 2: Loop extraction and sorting =================
for (cluster in target_clusters) {
  
  # --- A. Find Top Activity (already sorted by AUC descending) ---
  if (cluster %in% colnames(mean_auc_df)) {
    # Extract and maintain order (Rank 1 -> Rank N)
    top_activity <- rownames(mean_auc_df)[order(mean_auc_df[[cluster]], decreasing = TRUE)][1:top_n_cutoff]
  } else {
    top_activity <- c()
  }
  
  # --- B. Find Top Specificity (already sorted by RSS descending) ---
  if (cluster %in% colnames(rss_matrix)) {
    # Extract and maintain order
    top_rss <- rownames(rss_matrix)[order(rss_matrix[, cluster], decreasing = TRUE)][1:top_n_cutoff]
  } else {
    top_rss <- c()
  }
  
  # --- C. Take union and maintain "ranking" logic ---
  # unique() retains the order of first appearance.
  # So the result is: list high activity first (by AUC), then supplement with high specificity (by RSS)
  union_set <- unique(c(top_activity, top_rss))
  
  # Remove NAs (just in case)
  union_set <- union_set[!is.na(union_set)]
  
  # Store in List
  all_regulons_list[[cluster]] <- union_set
  
  cat("Extracted:", cluster, "- Count:", length(union_set), "\n")
}

# ================= STEP 3: Pad lengths and merge (Core step) =================

# 1. Find the number of Regulons in the longest subpopulation
max_length <- max(sapply(all_regulons_list, length))

# 2. Define a function: pad short vectors with empty strings "" until max_length is reached
pad_vector <- function(x, max_len) {
  c(x, rep("", max_len - length(x))) # Pad with empty strings at the end
}

# 3. Pad each subpopulation in the List and convert to Data Frame
# stringAsFactors = FALSE prevents characters from becoming factors
final_df <- data.frame(lapply(all_regulons_list, pad_vector, max_len = max_length), 
                       check.names = FALSE, 
                       stringsAsFactors = FALSE)

# 4. (Optional) Add a Rank column to easily see rankings
final_df <- cbind(Rank = 1:nrow(final_df), final_df)

# ================= STEP 4: Save results =================
output_filename <- "Cluster_Regulons_1231.csv"
write.csv(final_df, file = output_filename, row.names = FALSE, fileEncoding = "UTF-8")

cat("\n================================================\n")
cat("Success! CSV generated: ", output_filename, "\n")
cat("Format note: Each column is a subpopulation, sorted top to bottom by importance.\n")
#==========================================================


# ==========================================================
# ==========================================================
# Plot display graphs, scatter plots, important regulon CD4+ Treg (Th1-Like Transitional)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ==============================================================================
# 0. Define Morandi Palette
# ==============================================================================
morandi_red  <- "#A85751"  # Morandi Brick Red (for mechanistic genes)
morandi_blue <- "#5E8295"  # Morandi Slate Blue (for Top genes)
bg_color     <- "#333333"  # Deep Charcoal (for background points, softer than pure black)

# ==============================================================================
# 1. Prepare plotting data (follow your previous logic)
# ==============================================================================
target_name <- "CD4+ Treg (Th1-Like Transitional)"
if(!exists("rss_matrix") | !exists("mean_auc_df")) stop("Please run the calculation code first!")

common_genes <- intersect(rownames(rss_matrix), rownames(mean_auc_df))
df_scatter <- data.frame(
  Regulon = common_genes,
  RSS = rss_matrix[common_genes, target_name],
  AUC = mean_auc_df[common_genes, target_name],
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. Define highlights and labels (logic unchanged)
# ==============================================================================
# A. Auto-extract Top 5
top_rss_genes <- df_scatter %>% top_n(5, wt = RSS) %>% pull(Regulon)
top_auc_genes <- df_scatter %>% top_n(5, wt = AUC) %>% pull(Regulon)

# B. Manually specify mechanistic genes
mechanistic_genes <- c("HIF1A(+)","STAT1(+)", "NR3C1(+)")

# C. Group labeling
df_scatter$Label_Group <- "Other"
df_scatter$Label_Group[df_scatter$Regulon %in% c(top_rss_genes, top_auc_genes)] <- "Top Specific/Active"
valid_mechanistic <- intersect(mechanistic_genes, df_scatter$Regulon)
df_scatter$Label_Group[df_scatter$Regulon %in% valid_mechanistic] <- "Mechanistic Driver"

# D. Determine labels
genes_to_label <- unique(c(top_rss_genes, top_auc_genes, valid_mechanistic))

# ==============================================================================
# 3. Plotting (Color scheme updated)
# ==============================================================================
p_scatter <- ggplot(df_scatter, aes(x = AUC, y = RSS)) +
  
  # 1. Plot background points (darken color)
  # Use bg_color (#333333), set alpha to 0.4 to ensure overlaps are visible, but dark enough overall
  geom_point(data = subset(df_scatter, Label_Group == "Other"),
             color = bg_color, size = 1.2, alpha = 0.4) +
  
  # 2. Plot Top genes (Morandi Blue)
  geom_point(data = subset(df_scatter, Label_Group == "Top Specific/Active"),
             color = morandi_blue, size = 2.0, alpha = 0.9) +
  
  # 3. Plot Mechanistic genes (Morandi Red - slightly larger to stand out)
  geom_point(data = subset(df_scatter, Label_Group == "Mechanistic Driver"),
             color = morandi_red, size = 2.0) +
  
  # 4. Add text labels (font color follows point color)
  geom_text_repel(data = subset(df_scatter, Regulon %in% genes_to_label),
                  aes(label = Regulon, color = Label_Group), # Font color mapping
                  max.overlaps = Inf, 
                  box.padding = 0.6,
                  fontface = "bold",
                  bg.color = "white", # Add a white stroke to the font to prevent dark background points from interfering with reading
                  bg.r = 0.15,
                  show.legend = FALSE) +
  
  # 5. Apply Morandi color palette
  scale_color_manual(values = c("Mechanistic Driver" = morandi_red, 
                                "Top Specific/Active" = morandi_blue)) +
  
  theme_classic() +
  labs(title = paste0("Drivers of Fragility: ", target_name),
       subtitle = "Highlighting mTORC1-HIF1A Metabolic Reprogramming",
       x = "Mean Activity (AUC)", 
       y = "Specificity Score (RSS)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
        axis.title = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_blank()) # Remove legend title for simplicity

print(p_scatter)


# Show regulon differences between the two subpopulations using a Venn diagram
#install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

# ==============================================================================
# 1. Prepare data (keep unchanged)
# ==============================================================================
venn_data <- list(
  "Treg (Th1-Like)" = union_regulons_list[["CD4+ Treg (Th1-Like Transitional)"]],
  "T_ISG (STAT1+)"  = union_regulons_list[["CD4+ T_ISG (STAT1+)"]]
)

# ==============================================================================
# 2. Plotting (New: scale_x_continuous to expand canvas)
# ==============================================================================
p_venn <- ggVennDiagram(venn_data, 
                        label_alpha = 0,    
                        label_size = 4,     
                        set_size = 4) +     
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  
  # [Key fix]: Add 20% padding to both left and right to ensure long labels are not cut off
  scale_x_continuous(expand = expansion(mult = .2)) + 
  
  labs(title = "Regulon Overlap: Treg (Th1-Like) vs T_ISG",
       subtitle = "Shared IFN-response vs Unique Metabolic Drivers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

print(p_venn)

# ==============================================================================
# 3. (Optional) Print specific differentially expressed genes for use in your article
# ==============================================================================
unique_treg <- setdiff(venn_data[[1]], venn_data[[2]])
unique_isg  <- setdiff(venn_data[[2]], venn_data[[1]])
shared      <- intersect(venn_data[[1]], venn_data[[2]])

cat("\n=== Quick Check of Core Differences ===\n")
cat("🔴 Unique to Treg (Th1-Like) (Focus on HIF1A, REL, MYCN):\n")
print(unique_treg)
cat("\n🔵 Unique to T_ISG (Focus on NFE2L2, RUNX3):\n")
print(unique_isg)
cat("\n🟣 Shared by both (Focus on STAT1, IRF):\n")
print(shared)



# ==========================================================
# ==========================================================
# Plot display graph, scatter plot, important regulon CD4+ T_ISG (STAT1+)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ==============================================================================
# 1. Prepare data (specifically for T_ISG subpopulation)
# ==============================================================================
target_name <- "CD4+ T_ISG (STAT1+)" # Ensure name exactly matches

if(!exists("rss_matrix") | !exists("mean_auc_df")) stop("Please run the previous calculation code first!")

common_genes <- intersect(rownames(rss_matrix), rownames(mean_auc_df))
df_scatter_isg <- data.frame(
  Regulon = common_genes,
  RSS = rss_matrix[common_genes, target_name],
  AUC = mean_auc_df[common_genes, target_name],
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. Define highlight strategy (Focus: Antioxidant Resistance & Exhaustion)
# ==============================================================================
# A. Auto-extract Top 5
top_rss_genes <- df_scatter_isg %>% top_n(5, wt = RSS) %>% pull(Regulon)
top_auc_genes <- df_scatter_isg %>% top_n(5, wt = AUC) %>% pull(Regulon)

# B. [Key] Manually specify "Resistance/Exhaustion" mechanistic genes
# These genes explain why it is a Non-responder
mechanistic_genes <- c(
  "BCLAF1(+)",  # Core: Represents Pro-apoptotic Priming. This is a key molecule in the AICD execution phase, proving this group of cells is heading towards death
  "REL(+)",   # Core: High expression of REL proves the cells are in an "Activated but Dying" state.
  "STAT1(+)",   # Core: Chronic interferon signaling (interferon toxicity)
  "BACH1(+)",    # Feature: Differentiation Arrest. Explains why they have signals but no function. They are "locked" by these stemness/inhibitory TFs and cannot become real Effectors.
  "NFE2L2(+)"   # Feature: Represents Oxidative Stress. This is direct evidence of mitochondrial overload caused by the cell's inability to handle excessive IFNg stimulation.
)

# C. Group labeling (logic same as above)
df_scatter_isg$Label_Group <- "Other"
df_scatter_isg$Label_Group[df_scatter_isg$Regulon %in% c(top_rss_genes, top_auc_genes)] <- "Top Specific/Active"
valid_mechanistic <- intersect(mechanistic_genes, df_scatter_isg$Regulon)
df_scatter_isg$Label_Group[df_scatter_isg$Regulon %in% valid_mechanistic] <- "Mechanistic Driver"

# D. Determine labels
genes_to_label <- unique(c(top_rss_genes, top_auc_genes, valid_mechanistic))

# ==============================================================================
# 3. Plotting (Morandi color palette)
# ==============================================================================
# Define colors (slightly change the hue to distinguish from Treg)
# Mechanistic genes use "Morandi Rust Red" (#9E4838) -> Warn of resistance
# Top genes use "Morandi Deep Sea Blue" (#466986)
morandi_rust <- "#9E4838"
morandi_ocean <- "#466986"
bg_color <- "#333333"

p_scatter_isg <- ggplot(df_scatter_isg, aes(x = AUC, y = RSS)) +
  
  # Background points
  geom_point(data = subset(df_scatter_isg, Label_Group == "Other"),
             color = bg_color, size = 1.2, alpha = 0.4) +
  
  # Top genes
  geom_point(data = subset(df_scatter_isg, Label_Group == "Top Specific/Active"),
             color = morandi_ocean, size = 2.5, alpha = 0.9) +
  
  # Mechanistic genes (Focus)
  geom_point(data = subset(df_scatter_isg, Label_Group == "Mechanistic Driver"),
             color = morandi_rust, size = 3.5) +
  
  # Text labels
  geom_text_repel(data = subset(df_scatter_isg, Regulon %in% genes_to_label),
                  aes(label = Regulon, color = Label_Group),
                  max.overlaps = Inf, 
                  box.padding = 0.6,
                  fontface = "bold",
                  bg.color = "white", 
                  bg.r = 0.15,
                  show.legend = FALSE) +
  
  scale_color_manual(values = c("Mechanistic Driver" = morandi_rust, 
                                "Top Specific/Active" = morandi_ocean)) +
  
  theme_classic() +
  labs(title = paste0("Drivers of Resistance: ", target_name),
       subtitle = "Highlighting NRF2-Antioxidant & Blimp1-Exhaustion Axis",
       x = "Mean Activity (AUC)", 
       y = "Specificity Score (RSS)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
        legend.position = "bottom")

print(p_scatter_isg)

#================================Modified on 2026.1.29================================
# ==============================================================================
# 0. Define color scheme (Unified Red Style)
# ==============================================================================
morandi_red  <- "#A85751"  # Morandi Brick Red (for all core Regulons)
bg_color     <- "#333333"  # Deep Charcoal (for background points)

# ==============================================================================
# 1. Prepare plotting data
# ==============================================================================
target_name <- "CD4+ Treg (Th1-Like Transitional)"
if(!exists("rss_matrix") | !exists("mean_auc_df")) stop("Please run the calculation code first!")

common_genes <- intersect(rownames(rss_matrix), rownames(mean_auc_df))
df_scatter <- data.frame(
  Regulon = common_genes,
  RSS = rss_matrix[common_genes, target_name],
  AUC = mean_auc_df[common_genes, target_name],
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. Define highlights and labels (Merge lists, unified labeling)
# ==============================================================================

# Define all core genes to highlight (including Adaptation/Brake and Identity/Stemness)
# Exclude noise like POU1F1, HOXC11
highlight_genes <- c(
  # Adaptation & Brake (Mechanistic Core)
  "HIF1A(+)", "STAT1(+)", "NR3C1(+)", 
  # Identity & Stemness (Identity Maintenance)
  "LEF1(+)", "IKZF1(+)", "RUNX1(+)", "KLF12(+)", "FOXN3(+)"
)

# Initialize grouping
df_scatter$Label_Group <- "Other"

# Mark core genes
valid_highlight <- intersect(highlight_genes, df_scatter$Regulon)
df_scatter$Label_Group[df_scatter$Regulon %in% valid_highlight] <- "Key Regulon"

# ==============================================================================
# 3. Plotting (Unified Red)
# ==============================================================================
p_scatter <- ggplot(df_scatter, aes(x = AUC, y = RSS)) +
  
  # 1. Plot background points (Grey)
  geom_point(data = subset(df_scatter, Label_Group == "Other"),
             color = bg_color, size = 1.2, alpha = 0.3) +
  
  # 2. Plot Core genes (Unified Morandi Red)
  geom_point(data = subset(df_scatter, Label_Group == "Key Regulon"),
             color = morandi_red, size = 2.5) +  # Slightly increase point size
  
  # 3. Add text labels (Unified color, bold)
  geom_text_repel(data = subset(df_scatter, Regulon %in% valid_highlight),
                  aes(label = Regulon), 
                  color = morandi_red,      # Font color also uses red
                  max.overlaps = Inf, 
                  box.padding = 0.6,
                  point.padding = 0.4,
                  fontface = "bold",
                  bg.color = "white",       # White stroke, prevents background interference
                  bg.r = 0.15,
                  show.legend = FALSE) +
  
  theme_classic() +
  
  # 4. Titles and axis labels
  labs(title = "Drivers of Resilience: CD4+ Th1-like Treg",
       subtitle = "Coexistence of Adaptation (NR3C1/HIF1A) & Identity Maintenance",
       x = "Mean Activity (AUC)", 
       y = "Specificity Score (RSS)") +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none") # Legend not needed, since colors are unified

print(p_scatter)

#=========================TISG=================================================
# ==============================================================================
# 0. Define color scheme (Unified Red Style, echoing the Treg plot)
# ==============================================================================
morandi_red  <- "#A85751"  # Morandi Brick Red (for all core Regulons)
bg_color     <- "#333333"  # Deep Charcoal (for background points)

# ==============================================================================
# 1. Prepare data (For T_ISG subpopulation)
# ==============================================================================
target_name <- "CD4+ T_ISG (STAT1+)" # Ensure name exactly matches

if(!exists("rss_matrix") | !exists("mean_auc_df")) stop("Please run the calculation code first!")

common_genes <- intersect(rownames(rss_matrix), rownames(mean_auc_df))
df_scatter_isg <- data.frame(
  Regulon = common_genes,
  RSS = rss_matrix[common_genes, target_name],
  AUC = mean_auc_df[common_genes, target_name],
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. Define highlights (Only show core Regulons related to "Collapse/Drug Resistance")
# ==============================================================================

# Manually specify the 6 main "tragic" protagonists we discussed
highlight_genes <- c(
  "NR3C1(+)",   # [Core] Death Accelerator: Drives dysfunction synergistically with STAT1 in the absence of TGFb protection
  "BCLAF1(+)",  # [Outcome] AICD Execution: Pro-apoptotic priming, a sign of heading towards death
  "REL(+)",     # [Cause] Overactivation: Hyperactivation of NFkB pathway, a prerequisite for AICD
  "NFE2L2(+)",  # [Environment] Metabolic Collapse: Anti-oxidative stress (ROS), evidence of mitochondrial overload
  "STAT1(+)",   # [Source] IFN Imprint: Strong interferon stimulation signal
  "BACH1(+)"    # [State] Differentiation Arrest: Prevents cells from differentiating into effective effector cells
)

# Group labeling
df_scatter_isg$Label_Group <- "Other"
valid_highlight <- intersect(highlight_genes, df_scatter_isg$Regulon)
df_scatter_isg$Label_Group[df_scatter_isg$Regulon %in% valid_highlight] <- "Mechanistic Driver"

# ==============================================================================
# 3. Plotting (Unified Red Display)
# ==============================================================================
p_scatter_isg <- ggplot(df_scatter_isg, aes(x = AUC, y = RSS)) +
  
  # 1. Background points
  geom_point(data = subset(df_scatter_isg, Label_Group == "Other"),
             color = bg_color, size = 1.2, alpha = 0.3) +
  
  # 2. Mechanistic genes (Morandi Red)
  geom_point(data = subset(df_scatter_isg, Label_Group == "Mechanistic Driver"),
             color = morandi_red, size = 2.5) +
  
  # 3. Text labels (Red bold)
  geom_text_repel(data = subset(df_scatter_isg, Regulon %in% valid_highlight),
                  aes(label = Regulon),
                  color = morandi_red,      # Font uniformly uses red
                  max.overlaps = Inf, 
                  box.padding = 0.6,
                  fontface = "bold",
                  bg.color = "white",       # White stroke
                  bg.r = 0.15,
                  show.legend = FALSE) +
  
  theme_classic() +
  
  # 4. Title update: Emphasize AICD and Collapse
  labs(title = paste0("Drivers of Resistance: ", target_name),
       subtitle = "Hyperactivation (REL), Stress (NR3C1/NFE2L2) & AICD (BCLAF1)",
       x = "Mean Activity (AUC)", 
       y = "Specificity Score (RSS)") +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none") # Remove legend

print(p_scatter_isg)

















# View what the target genes of these upregulated Regulons are
library(data.table)
library(stringr)

# 1. Read reg.csv
reg_data <- fread("reg.csv")

# 2. Define a function to parse target genes
# The TargetGenes column in reg.csv looks like this: [('GeneA', 0.5), ('GeneB', 0.8)...]
get_targets <- function(tf_name, reg_df) {
  # Find the row for this TF (there might be multiple rows for different Motifs, usually choose the one with the highest Enrichment)
  tf_row <- reg_df[reg_df$TF == tf_name, ]
  
  if (nrow(tf_row) == 0) return(NULL)
  
  # Take the first row (usually the most significant)
  targets_str <- tf_row$TargetGenes[1]
  
  # Clean string, extract gene names
  # Remove parentheses, quotes, fractions, leave only gene names
  # This is a simple regex extraction method
  genes <- str_extract_all(targets_str, "'[^']+'")[[1]]
  genes <- gsub("'", "", genes) # Remove single quotes
  
  return(genes)
}

# --- Practical operation ---

# Suppose you see "STAT1" is the Top TF of T_ISG in the plot from step 3
target_genes_stat1 <- get_targets("STAT1", reg_data)

# Print the top 20 target genes
print(head(target_genes_stat1, 20))

# Advanced: Do a GO enrichment analysis to see what these target genes do
library(clusterProfiler)
library(org.Hs.eg.db) # Assuming human

ego <- enrichGO(gene = target_genes_stat1,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")
barplot(ego, showCategory = 10, title = "STAT1 Targets GO Enrichment")