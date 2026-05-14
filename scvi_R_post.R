## Change system error messages to English
Sys.setenv(LANGUAGE = "en")
## Disable conversion to factors
options(stringsAsFactors = FALSE)
## Clear environment
rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)
getwd()
setwd("D:/Data_us/data/res/")
#===============================================================================
### Step 1: Read data exported from Python
# 1. Set file path (Please modify to the actual folder path where you saved the CSV)
# If you used D:/Data_us/data/combined_data/scvi_export_for_r in Python earlier
scvi_path <- "D:/Data_us/data/combined_data/scvi_export_for_r/scvi_embeddings.csv"

# 2. Read scVI coordinates file
# row.names = 1 is very important because the first column is the cell name
scvi_embeddings <- read.csv(scvi_path, row.names = 1)

# 3. Check if it's read correctly
print(dim(scvi_embeddings)) # Should show (around 210,000, 10)
head(scvi_embeddings)

#===============================================================================
# Step 2: Import data into Seurat object (Revised version - Auto-align cells)

# 1. Restore object and clear memory
# Whatever state your sce.combined is in now, let's reset it first
sce.combined <- sce_combined_backup # No matter what object it is, assign it and then delete
rm(sce_combined_backup)

# 2. [Key step] Intersect and filter Seurat object
# Some low-quality cells were filtered during Python processing, they must be synchronously removed in the R object
common_cells <- intersect(rownames(scvi_embeddings), Cells(sce.combined))

print(paste0("Original Seurat cell count: ", length(Cells(sce.combined))))
print(paste0("scVI result contains cell count: ", length(common_cells)))

# Slim down the Seurat object, keeping only the cells present in the scVI results
sce.combined <- subset(sce.combined, cells = common_cells)

# 3. Ensure the order of scVI data matches the slimmed-down Seurat object exactly
scvi_embeddings <- scvi_embeddings[Cells(sce.combined), ]

# 4. Check matching again (Double insurance)
if(all(rownames(scvi_embeddings) == Cells(sce.combined))) {
  message("✅ Cell names matched perfectly, low-quality cells synchronously removed!")
} else {
  stop("❌ Cell names still do not match, please check code logic.")
}

# 5. Put scVI coordinates into Seurat
# Seurat requires column names to match the key (e.g., if key="scVI_", column names must be scVI_1, scVI_2...)
colnames(scvi_embeddings) <- paste0("scVI_", 1:ncol(scvi_embeddings))
sce.combined[["scvi"]] <- CreateDimReducObject(
  embeddings = as.matrix(scvi_embeddings), 
  key = "scVI_", 
  assay = DefaultAssay(sce.combined)
)

print("scVI coordinates successfully injected!")

# ==============================================================================
# Step 3: Downstream analysis (Clustering and UMAP)


print("Performing UMAP dimensionality reduction (based on scVI)...")
# Note: RunPCA and ScaleData are not needed, use scVI directly
sce.combined <- RunUMAP(sce.combined, reduction = "scvi", dims = 1:10)

print("Finding neighbors (FindNeighbors)...")
sce.combined <- FindNeighbors(sce.combined, reduction = "scvi", dims = 1:10)

print("Clustering (FindClusters)...")
# Resolution is adjustable: 0.5 is default, try 0.8 or 1.0 for large cell numbers
sce.combined <- FindClusters(sce.combined, resolution = 0.8)

print("✅ All analyses completed! You can now use DimPlot(sce.combined, reduction = 'umap') to view the results.")

# Plotting
# Figure 1: View cluster distribution (with numeric labels)
p1 <- DimPlot(sce.combined, reduction = "umap", label = TRUE, label.size = 4, raster = FALSE) + 
  ggtitle(paste("scVI Clustering (Res 0.8) - 36 Clusters")) +
  theme(plot.title = element_text(hjust = 0.5))

# Figure 2: View batch integration effect (colored by sample)
# If 27 samples are perfectly mixed together, batch effect removal is perfect
p2 <- DimPlot(sce.combined, reduction = "umap", group.by = "orig.ident", raster = FALSE) + 
  ggtitle("Batch Integration Check") + 
  theme(plot.title = element_text(hjust = 0.5)) # Too many samples, hide legend first to avoid blocking the plot

# Display plots
p1
p2
# ===============================================================================
# Downstream annotation
library(dplyr)
message("Starting to calculate Marker genes for all subclusters...")
message("To speed up, a maximum of 500 cells per cluster are sampled for calculation.")

# Ensure default assay is RNA counts, so differential expression calculation is most accurate
DefaultAssay(sce.combined) <- "RNA"

# Find Markers for all clusters
# min.pct = 0.25: Gene expressed in at least 25% of cells in the potential target cluster
# logfc.threshold = 0.25: Log fold change must be at least 0.25
# only.pos = TRUE: Only find highly expressed Markers
markers_all <- FindAllMarkers(
  sce.combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500, # [Key acceleration parameter]
  test.use = "wilcox",       # Default test method used
  verbose = TRUE
)

message("Calculation completed! Organizing results...")

# Filter the top 10 most significant Markers for each cluster (sorted by adjusted p-value)
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Strongly recommend saving this result so you don't have to calculate it next time
setwd("D:/Data_us/data/combined_data/scvi_export_for_r/")
write.csv(markers_all, "all_cluster_markers.csv")
write.csv(top10_markers, "top10_cluster_markers_summary.csv")

message("✅ Marker gene table has been saved as a CSV file.")
# You can use View(top10_markers) in RStudio for a preliminary look

# ===============================================================================
# Define Marker gene list for common cell types
# Note: Gene names must be all uppercase because they are usually uppercase in your scRNA data
marker_list <- list(
  # --- 1. Tumor cells (Melanoma specific) ---
  # MITF: Key transcription factor for melanocyte development
  # PMEL (gp100), MLANA (MART-1), TYR (Tyrosinase): Melanin synthesis pathway genes, classic clinical Markers
  "Melanoma_Tumor" = c("MITF", "PMEL", "MLANA", "TYR", "S100B"),
  
  # --- 2. Epithelial/Keratinocytes (If samples contain normal skin tissue) ---
  # KRT14, KRT5: Basal cell layer Markers
  # KRT10, KRT1: Spinous layer/Differentiated cell layer Markers
  "Epithelial/Keratinocytes" = c("KRT14", "KRT5", "KRT10", "KRT1"),
  
  # --- 3. Stromal Cells ---
  # Fibroblasts: COL1A1 (Collagen), DCN, PDGFRB
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"),
  # Endothelial cells: PECAM1 (CD31), VWF, CLDN5
  "Endothelial" = c("PECAM1", "VWF", "CLDN5", "FLT1"),
  # Pericytes/Smooth Muscle Cells (SMC): ACTA2 (SMA), RGS5
  "Pericytes/SMC" = c("ACTA2", "RGS5", "TAGLN"),
  
  # --- 4. Lymphoid immune cells (Lymphoid) ---
  # T Cells (Total T cells): CD3D, CD3E
  "T_Cells_General" = c("CD3D", "CD3E"),
  # CD8+ T Cells (Cytotoxic T cells): CD8A, CD8B, GZMB (Granzyme B, activation marker)
  "CD8_T_Cytotoxic" = c("CD8A", "CD8B", "GZMB", "NKG7"),
  # CD4+ T Cells (Helper T cells): CD4, IL7R
  "CD4_T_Helper" = c("CD4", "IL7R"),
  # Tregs (Regulatory T cells): FOXP3, IL2RA (CD25) - Suppress anti-tumor immunity
  "Tregs" = c("FOXP3", "IL2RA", "CTLA4"),
  # NK Cells (Natural Killer cells): GNLY, KLRD1 (NKG7 and CD8 are also high in NK, need comprehensive judgment)
  "NK_Cells" = c("GNLY", "KLRD1", "NCAM1"),
  # B Cells: MS4A1 (CD20), CD79A
  "B_Cells" = c("MS4A1", "CD79A"),
  # Plasma Cells: JCHAIN, MZB1, IGHG1 (Produce antibodies)
  "Plasma_Cells" = c("JCHAIN", "MZB1", "IGHG1"),
  
  # --- 5. Myeloid immune cells (Myeloid) ---
  # Macrophages/Monocytes: CD68, CD163 (M2-biased anti-inflammatory), CSF1R, LYZ
  "Myeloid_Macro_Mono" = c("CD68", "CD163", "CSF1R", "LYZ", "AIF1"),
  # Dendritic Cells: CD1C (cDC2), LILRA4 (pDC), CLEC9A (cDC1)
  "Dendritic_Cells" = c("CD1C", "LILRA4", "CLEC9A", "HLA-DRA"),
  # Mast Cells: TPSAB1 (Tryptase), CPA3
  "Mast_Cells" = c("TPSAB1", "CPA3"),
  
  # --- 6. Proliferating Cells (Non-cell-type specific, indicating active division) ---
  # MKI67 (Ki-67), TOP2A
  "Proliferation_Cycle" = c("MKI67", "TOP2A")
)

message("Marker list preparation completed.")

# ===============================================================================
# 1. Organize Marker list
# Flatten the list into a vector and remove duplicates
genes_to_plot <- unique(unlist(marker_list))

# 2. [Important] Check if these genes exist in your data
# Some genes might have been filtered out during preprocessing due to extremely low expression
genes_present <- intersect(genes_to_plot, rownames(sce.combined))
genes_missing <- setdiff(genes_to_plot, rownames(sce.combined))

if(length(genes_missing) > 0) {
  message("⚠️ Note: The following genes were not found in your data and will not be displayed in the plot:")
  print(genes_missing)
}

# 3. Plot DotPlot
# group.by = "seurat_clusters" ensures grouping by clustering results
# assay = "RNA" ensures using raw counts to reflect true expression abundance
dot_plot <- DotPlot(sce.combined, 
                    features = genes_present, 
                    group.by = "seurat_clusters",
                    assay = "RNA", 
                    cols = c("lightgrey", "red"), # Color gradient: light grey to red
                    cluster.idents = FALSE) +     # Maintain cluster numbering order (0, 1, 2...)
  RotatedAxis() + # Rotate X-axis gene names to prevent overlap
  theme(axis.text.x = element_text(size = 9), # Adjust font size
        axis.text.y = element_text(size = 10)) +
  labs(title = "Marker Expression across 36 Clusters",
       x = "Markers", y = "Cluster ID")

# 4. Display and save the plot
# The plot will be very long due to 36 clusters and dozens of genes, recommend saving as a large PDF or PNG for viewing
print(dot_plot)

ggsave("DotPlot_CellType_Markers.pdf", plot = dot_plot, width = 18, height = 12)
ggsave("DotPlot_CellType_Markers.png", plot = dot_plot, width = 18, height = 12)
message("✅ DotPlot has been saved, please open the PDF/PNG file for detailed review.")

# ===============================================================================
# Annotate based on DotPlot results
# 1. Create a mapping vector from ID to cell type
# Revised cell type annotation (Revised Annotation)

new.cluster.ids <- c(
  # --- Tumor cells (Melanoma) ---
  # Features: MITF+, PMEL+, MLANA+
  "0" = "Melanoma",
  "1" = "Melanoma",
  "2" = "Melanoma",
  "3" = "Melanoma",
  "5" = "Melanoma",
  "7" = "Melanoma",
  "9" = "Melanoma",
  "18" = "Melanoma",
  "19" = "Melanoma",
  "25" = "Melanoma (MHC-II+)", # HLA-DRA high
  
  # --- T Cells & NK Cells ---
  "4" = "CD4+ T Cells (Naive/Mem)", # IL7R+, CD4+, CTLA4-
  "10" = "Tregs",                   # CTLA4+, CD4+ (Key correction)
  "6" = "CD8+ T Cells",             # CD8A+, NKG7+
  "29" = "NK Cells (GNLY+)",        # GNLY+, NKG7+, CD3-
  "28" = "NK Cells (GZMB+)",        # GZMB+, PECAM1- (Key correction: Previously misidentified as endothelial)
  
  # --- B Cells & Plasma Cells ---
  "8" = "B Cells",                  # MS4A1+
  "24" = "B Cells",
  "14" = "Plasma Cells",            # JCHAIN+, MZB1+
  
  # --- Myeloid Cells ---
  "16" = "Macrophages (M2-like)",   # CD163+, CD68+
  "32" = "Monocytes",               # LYZ+, HLA-DRA+, CD163- (Key correction: Monocytes)
  "17" = "Dendritic Cells (cDC2)",  # CD1C+, HLA-DRA+
  "33" = "Mast Cells",              # TPSAB1+, CPA3+
  
  # --- Stroma and Vasculature ---
  "11" = "Fibroblasts",             # COL1A1+
  "12" = "Fibroblasts",
  "15" = "Fibroblasts (Matrix)",    # DCN++, LUM++
  "26" = "Fibroblasts",
  "35" = "Pericytes/SMC",           # ACTA2+, TAGLN+ (Key correction: Previously misidentified as Unknown)
  "13" = "Endothelial",             # PECAM1+, VWF+
  "27" = "Endothelial",
  
  # --- Epithelial/Keratinocytes ---
  "21" = "Keratinocytes (Basal)",   # KRT14+
  "23" = "Keratinocytes (Cycling)", # KRT14+, MKI67+
  "30" = "Keratinocytes (Diff)",    # KRT10+
  "31" = "Keratinocytes (Diff)",
  
  # --- Other/Doublets/Low Quality ---
  "20" = "Proliferating Cells",     # MKI67+, TOP2A+ (Could be proliferating tumor or immune cells)
  "22" = "Doublets (Epi+Plasma)",   # High KRT+ & JCHAIN+ simultaneously
  "34" = "Unknown"                  # Signal too weak, no distinct Markers
)

# Execute renaming
Idents(sce.combined) <- "RNA_snn_res.0.8"
sce.combined <- RenameIdents(sce.combined, new.cluster.ids)
sce.combined$cell_type <- Idents(sce.combined)

# Plot again
DimPlot(sce.combined, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle("Revised Cell Type Annotations")


# ===============================================================================
# Found many doublets, remove doublets
library(scDblFinder)
library(SingleCellExperiment)

message("Converting Seurat object to SingleCellExperiment format...")
# Ensure using raw RNA count data for calculation
DefaultAssay(sce.combined) <- "RNA"

# 1. Convert to SCE object
sce <- as.SingleCellExperiment(sce.combined)

# 2. Run scDblFinder
# [Key adaptation] Added samples="orig.ident" here, consistent with your original function logic
# (Originally run on single sample, now run separately on each single sample in the merged object)
message("Running scDblFinder (this may take a few minutes)...")
sce <- scDblFinder(sce, samples = "orig.ident") 

# 3. Extract results
message("Extracting doublet scores...")
# Extract doublet results from SCE (refer to your function logic)
df_result <- as.data.frame(colData(sce)[, c("scDblFinder.class", "scDblFinder.score")])

# 4. Ensure order matching (Seurat and SCE conversions usually keep order consistent, but adding match is safer)
df_result <- df_result[match(rownames(sce.combined@meta.data), rownames(df_result)), ]

# 5. Add to Seurat v5 object's meta.data
sce.combined$doublet_class <- df_result$scDblFinder.class
sce.combined$doublet_score <- df_result$scDblFinder.score

# Count how many doublets were found
table(sce.combined$doublet_class)

# Set save path (refer to your function structure, we create a folder)
if(!dir.exists("qcplot")) dir.create("qcplot")

message("Plotting doublet distribution...")

# Plot UMAP: doublets marked in red
pdf("./qcplot/UMAP_pre_doublet_removal_Combined.pdf", width = 10, height = 8)
print(DimPlot(sce.combined, reduction = "umap", group.by = "doublet_class", cols = c("#D51F26", "#272E6A")) + 
        ggtitle("scDblFinder: Singlet vs Doublet"))
dev.off()
getwd()

# Plot another image to compare our previous Cluster 22 and Cluster 20
# See if the doublets identified by scDblFinder are located in these clusters
p_compare <- DimPlot(sce.combined, reduction = "umap", group.by = "doublet_class", pt.size = 0.5)
p_compare

message("Removing doublets...")

# 1. Filter Singlets (refer to your function: subset = doublet_class != "doublet")
sce_clean <- subset(sce.combined, subset = doublet_class == "singlet" )

# 2. Plot UMAP after removal
pdf("./qcplot/UMAP_post_doublet_removal_Combined.pdf", width = 10, height = 8)
print(DimPlot(sce_clean, reduction = "umap")) + ggtitle("Cleaned Data (Singlets Only)")
dev.off()

# 3. Use DietSeurat to clean data (fully refer to your function)
message("Performing DietSeurat cleaning...")
sce_clean <- DietSeurat(sce_clean, 
                        assays = "RNA",                  # Keep only RNA assay
                        layers = c("counts", "data"),    # Keep counts and data
                        dimreducs = c("pca", "umap", "scvi"), # Keep dimensional reductions we will use (note I added scvi)
                        features = NULL                  # Keep all genes
)


# 4. Update final object
sce.combined <- sce_clean
rm(sce_clean) # Free memory
gc()# Toss to Gemini for a glance at the machine's running state and condition.

# Rerun UMAP (only for Singlets)
sce.combined <- RunUMAP(sce.combined, reduction = "scvi", dims = 1:10)

# Plot and check (Cluster 22 should have completely disappeared)
DimPlot(sce.combined, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Final Clean UMAP (Post-Doublet Removal)")

sce.combined <- FindNeighbors(sce.combined, reduction = "scvi", dims = 1:10)
# Resolution is adjustable: 0.5 is default, try 0.8 or 1.0 for large cell numbers
sce.combined <- FindClusters(sce.combined, resolution = 0.5)

rm(markers_all,marker_list)
rm(top10_markers)

# Figure 1: View cluster distribution (with numeric labels)
p3 <- DimPlot(sce.combined, reduction = "umap", label = TRUE, label.size = 4, raster = FALSE) + 
  ggtitle(paste("scVI Clustering (Res 0.5) - 28 Clusters")) +
  theme(plot.title = element_text(hjust = 0.5))
p3

# ===============================================================================
# Re-annotate
# Ensure default assay is RNA counts, so differential expression calculation is most accurate
sce.combined$seurat_clusters <- Idents(sce.combined)
DefaultAssay(sce.combined) <- "RNA"
markers_all <- FindAllMarkers(
  sce.combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500, # [Key acceleration parameter]
  test.use = "wilcox",       # Default test method used
  verbose = TRUE
)

# Filter the top 10 most significant Markers for each cluster (sorted by adjusted p-value)
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# group.by = "seurat_clusters" ensures grouping by clustering results
# assay = "RNA" ensures using raw counts to reflect true expression abundance
dot_plot <- DotPlot(sce.combined, 
                    features = genes_present, 
                    group.by = "seurat_clusters",
                    assay = "RNA", 
                    cols = c("lightgrey", "red"), # Color gradient: light grey to red
                    cluster.idents = FALSE) +     # Maintain cluster numbering order (0, 1, 2...)
  RotatedAxis() + # Rotate X-axis gene names to prevent overlap
  theme(axis.text.x = element_text(size = 9), # Adjust font size
        axis.text.y = element_text(size = 10)) +
  labs(title = "Marker Expression across 28 Clusters",
       x = "Markers", y = "Cluster ID")
dot_plot
p4 <- DimPlot(sce.combined, reduction = "umap", 
              group.by = "orig.ident",   # Core modification: Specify grouping/coloring by orig.ident
              label = FALSE,             # Recommend changing to FALSE: Samples are usually seen in the legend, labeling on the plot gets messy
              raster = FALSE,
              shuffle = TRUE) +          # Recommend adding: Shuffle the drawing order of points to avoid one sample blocking others
  ggtitle("UMAP colored by Sample (orig.ident)") +
  theme(plot.title = element_text(hjust = 0.5))

p4


# ==============================================================================
# Final version: Custom annotation and filtering

# 1. Prepare annotation vector (0-27 full coverage)
new.cluster.ids <- c(
  "0" = "Melanoma",
  "1" = "T Cells",
  "2" = "Melanoma",
  "3" = "Melanoma (HLA-DRA+)",
  "4" = "Myeloid",
  "5" = "B Cells",
  "6" = "Melanoma",
  "7" = "T Cells",
  "8" = "Melanoma",
  "9" = "Melanoma",
  "10" = "Endothelial",
  "11" = "Fibroblasts",
  "12" = "Plasma Cells",
  "13" = "Fibroblasts",
  "14" = "Myeloid",
  "15" = "Myeloid",
  "16" = "Melanoma",
  "17" = "Low Quality",           # To be deleted
  "18" = "Epithelial/Keratinocytes",
  "19" = "Pericytes",
  "20" = "B Cells",
  "21" = "T Cells",
  "22" = "Myeloid",              # [Corrected]
  "23" = "Low Quality",           # To be deleted
  "24" = "Endothelial",
  "25" = "Mast Cells",
  "26" = "Low Quality",           # [Corrected] Marked to be deleted
  "27" = "Epithelial/Keratinocytes"
)

# 2. Perform renaming
Idents(sce.combined) <- "seurat_clusters"

# This step converts Cluster IDs into your biological names
sce.combined <- RenameIdents(sce.combined, new.cluster.ids)

# 3. Store annotations into the cell_type column of metadata
sce.combined$cell_type <- Idents(sce.combined)

# Filter unneeded subclusters (Filter 17, 23, 26)
message(paste0("Cell count before filtering: ", ncol(sce.combined)))

# Delete all cells marked as "Low Quality" (i.e., original 17, 23, 26)
sce.final <- subset(sce.combined, idents = "Low Quality", invert = TRUE)

message(paste0("Cell count after filtering: ", ncol(sce.final)))

# Plot and save

# Plot final UMAP
p_final <- DimPlot(sce.final, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle("Final Manual Annotation (Cleaned)") +
  theme(legend.position = "right")

print(p_final)

# Count the number of each type of cell
print(table(Idents(sce.final)))

# Save final results
saveRDS(sce.final, "sce_final_manual_annotated.rds")
message("✅ All done! Cluster 22 has been grouped into NK, clusters 26, 17, 23 have been deleted.")
getwd()
table(sce.final$orig.ident)