library(Seurat)
library(ggplot2)
library(patchwork) # For patchwork/layout

# 1. Put your objects into a list and name them (for convenient subsequent labeling)
obj_list <- list(
  "Sample_1" = test_data,
  "Sample_T5" = test_data_T5,
  "Sample_T13" = test_data_T13,
  "Sample_T14" = test_data_T14
)

# =======================================================
# Mode A: View H&E image only (hide sequencing spots)
# =======================================================
# alpha = 0 makes sequencing spots fully transparent, stroke = 0 removes point borders
plots_he_only <- lapply(names(obj_list), function(x) {
  SpatialDimPlot(obj_list[[x]], alpha = 0, stroke = 0) + 
    ggtitle(paste0(x, " - H&E Only")) + 
    theme(legend.position = "none") # Legend not needed for viewing image only
})

# Patchwork display
wrap_plots(plots_he_only, ncol = 2)

# =======================================================
# Mode B: View H&E image + clustering spots (Routine check)
# =======================================================
# alpha controls point transparency (0.1-1), pt.size.factor controls point size
plots_clusters <- lapply(names(obj_list), function(x) {
  SpatialDimPlot(obj_list[[x]], 
                 alpha = 0.6,          # Points slightly transparent to reveal background
                 pt.size.factor = 1.6) + # Adjust point size as needed
    ggtitle(paste0(x, " - Clusters"))
})

# Patchwork display
wrap_plots(plots_clusters, ncol = 2)

library(Seurat)
library(SeuratDisk)

# =======================================================
# 1. Export single-cell reference
# =======================================================
# Assuming your single-cell object is named seurat_obj (containing CD8 Trm and Th1-like Treg annotations)
# Ensure your annotation column is named "cell_type" (or another name you specify, please note it down)
seurat_obj$cell_type <- seurat_obj$detailed_type_NichNet1 # For example, using your latest annotation column

# Convert to h5Seurat
SaveH5Seurat(seurat_obj, filename = "sc_reference.h5Seurat", overwrite = TRUE)
# Convert to h5ad
Convert("sc_reference.h5Seurat", dest = "h5ad", overwrite = TRUE)

# =======================================================
# 2. Export spatial data (Visium)
# =======================================================
library(Seurat)
library(SeuratDisk)
library(png)

ExportWithImage_Robust <- function(obj, sample_name) {
  
  print(paste("=== Processing sample:", sample_name, "==="))
  
  # -------------------------
  # 1. Ensure image can be found
  # -------------------------
  # Automatically look for the first image slot
  if (length(obj@images) == 0) stop("Error: No images found in this object!")
  img_key <- names(obj@images)[1]
  print(paste("  - Detected image Key:", img_key))
  
  spatial_obj <- obj@images[[img_key]]
  
  # -------------------------
  # 2. Save PNG image
  # -------------------------
  img_matrix <- spatial_obj@image
  if (is.null(img_matrix)) stop("Error: Image matrix is empty")
  
  png_filename <- paste0(sample_name, "_image.png")
  print(paste("  - Saving image:", png_filename))
  
  # Simple normalization to prevent errors
  if (max(img_matrix) > 1) {
    img_matrix <- img_matrix / 255
  }
  png::writePNG(img_matrix, target = png_filename)
  
  # -------------------------
  # 3. Save scale factors
  # -------------------------
  scale_factors <- spatial_obj@scale.factors
  sf_df <- data.frame(
    spot_diameter_fullres = scale_factors$spot,
    fiducial_diameter_fullres = scale_factors$fiducial,
    hires_scalef = scale_factors$hires,
    lowres_scalef = scale_factors$lowres
  )
  write.csv(sf_df, file = paste0(sample_name, "_scalefactors.csv"), row.names = FALSE)
  
  # -------------------------
  # 4. Extract coordinates (the most error-prone step, we use the safest method)
  # -------------------------
  print("  - Extracting coordinates...")
  
  # Extract coordinates using the specified image key
  coords <- GetTissueCoordinates(obj, image = img_key)
  
  # [Key fix] Do not use AddMetaData, manually assign directly instead
  # And automatically detect column names (Seurat V4 is imagerow/imagecol, V5 might be x/y)
  
  if ("imagerow" %in% colnames(coords)) {
    # Seurat V4 style: row=y, col=x
    obj$spatial_y <- coords[, "imagerow"]
    obj$spatial_x <- coords[, "imagecol"]
  } else if ("y" %in% colnames(coords)) {
    # Seurat V5 style
    obj$spatial_y <- coords[, "y"]
    obj$spatial_x <- coords[, "x"]
  } else {
    # Fallback plan: blindly take the first two columns (usually 1st column is y/row, 2nd is x/col, or vice versa)
    # Visium order here doesn't affect running Cell2Location, only plotting direction, can be transposed in Python
    print("  - Warning: Standard column names not recognized, using columns 1 and 2 as coordinates")
    obj$spatial_y <- coords[, 1]
    obj$spatial_x <- coords[, 2]
  }
  
  print("  - Coordinates extracted successfully, preparing to diet...")
  
  # -------------------------
  # 5. Diet and export h5ad
  # -------------------------
  # Clean up list columns in metadata (to prevent HDF5 errors)
  for (col in colnames(obj@meta.data)) {
    if (is.list(obj@meta.data[[col]])) obj@meta.data[[col]] <- NULL
  }
  
  # Use DietSeurat and [forcefully] remove images to prevent CreateGroup errors
  obj_diet <- DietSeurat(obj, assays = "Spatial", dimreducs = NULL, graphs = NULL)
  obj_diet@images <- list() 
  
  # Compatible with V5 Assay
  if (inherits(obj_diet[["Spatial"]], "Assay5")) {
    obj_diet[["Spatial"]] <- as(obj_diet[["Spatial"]], "Assay")
  }
  
  h5file <- paste0(sample_name, ".h5Seurat")
  print("  - Writing h5Seurat...")
  SaveH5Seurat(obj_diet, filename = h5file, overwrite = TRUE)
  
  print("  - Converting to h5ad...")
  Convert(h5file, dest = "h5ad", overwrite = TRUE)
  
  print(paste("=== Sample processing complete!", sample_name, "==="))
}

# ================= Rerun =================
ExportWithImage_Robust(test_data_T5, "spatial_T5")
ExportWithImage_Robust(test_data_T14, "spatial_T14")

# =======================================================
# Downsample all cells to obtain h5ad file
library(Seurat)
library(SeuratDisk)
library(dplyr)

# ==============================================================================
# 1. Prepare data and define paths
# ==============================================================================
# Assuming your large object is called sce_11_downstream
# Assuming your fine T cell object is called t_cells_clean_final_annotated3

# !!! Replace original_annotation here with the column in your sce_11_downstream
# storing rough annotations like "Tumor", "Fibroblast", "B cell" !!!
coarse_label_col <- "cell_type"  # Please modify according to actual situation, e.g., "seurat_clusters" or "major_type"

# Set downsampling threshold
sample_size <- 500 

# ==============================================================================
# 2. Merge annotation information (via Barcode matching)
# ==============================================================================
message("Start merging annotations...")

# Step 1: Initialize a new column and fill in the original rough annotations first
# [Core modification]: Must force convert to character using as.character to prevent Factor errors
sce_11_downstream$C2L_final_annotation <- as.character(sce_11_downstream@meta.data[[coarse_label_col]])

# Extract fine annotations of T cells
# Ensure the rownames (barcodes) of both objects are consistent
t_cell_barcodes <- intersect(colnames(sce_11_downstream), colnames(t_cells_clean_final_annotated3))

if(length(t_cell_barcodes) == 0) {
  stop("Error: Cell Barcodes of the two objects do not overlap! Please check if the Barcode prefixes have been changed.")
}

# Extract fine annotation vector and ensure it is also character
detailed_labels <- as.character(t_cells_clean_final_annotated3@meta.data[t_cell_barcodes, "detailed_type_NichNet1"])
names(detailed_labels) <- t_cell_barcodes

# Overwrite T cell annotations
sce_11_downstream$C2L_final_annotation[names(detailed_labels)] <- detailed_labels

# ==============================================================================
# 2. Clean and check
# ==============================================================================
# Set Idents to the new annotation
Idents(sce_11_downstream) <- "C2L_final_annotation"

message("Cell type statistics after merging:")
print(table(Idents(sce_11_downstream)))

# [Optional] If there is still a major class called "T cells" at this point (meaning some T cells are not subdivided),
# To avoid affecting Cell2location training, it is recommended to delete it or classify it as "Other T"
# Only run the line below if you see "T cells" remaining in the table output:
# sce_11_downstream <- subset(sce_11_downstream, idents = "T cells", invert = TRUE)

# ==============================================================================
# 3. Downsampling - Target: 500/class
# ==============================================================================
target_cells <- 500
message(paste0("Performing stratified downsampling, keeping up to per class: ", target_cells))

# Downsample based on current Idents (i.e., C2L_final_annotation)
sce_ref_final <- subset(sce_11_downstream, downsample = target_cells)

message("Number of cells finally used for Reference after downsampling:")
print(table(Idents(sce_ref_final)))

# Remove impure cells
# 1. Set Idents (ensure it's correct)
Idents(sce_ref_final) <- "C2L_final_annotation"

# 2. Remove "T Cells" (broad label) and "CD4+ Tem" (too few < 50)
# Note: Although Mast Cells (77) and Pericytes (84) are few, because they are vastly different from T cells/tumors,
# their features are distinct, so it is usually fine to keep them. But 17 is really too few.
cells_to_remove <- c("T Cells", "CD4+ Tem")

message("Removing interfering categories: ", paste(cells_to_remove, collapse = ", "))
sce_ref_clean <- subset(sce_ref_final, idents = cells_to_remove, invert = TRUE)

# 3. Check again
message("Final list used for Reference:")
print(table(Idents(sce_ref_clean)))

# ==============================================================================
# 4. Generate h5ad file
# ==============================================================================
# Use DietSeurat to diet, keep only Counts (Raw data is the core of C2L)
sce_export <- DietSeurat(
  sce_ref_final, 
  counts = TRUE, 
  data = TRUE, 
  scale.data = FALSE,
  dimreducs = NULL, 
  graphs = NULL
)

# Path settings (avoid paths with Chinese characters)
save_path <- "D:/Bioinfo_Temp/reference_melanoma_C2L.h5seurat"
h5ad_output <- "D:/Bioinfo_Temp/reference_melanoma_C2L.h5ad"

# Ensure directory exists
if(!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)

message("Saving h5seurat...")
SaveH5Seurat(sce_export, filename = save_path, overwrite = TRUE)

message("Converting to h5ad...")
Convert(save_path, dest = "h5ad", overwrite = TRUE)

message("All done! Output file is located at: ", h5ad_output)