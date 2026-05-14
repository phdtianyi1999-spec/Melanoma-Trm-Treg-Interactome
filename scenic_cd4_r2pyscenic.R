library(SCopeLoomR)
library(Seurat)

# 1. Extract the matrix
# This step will pop up an "allocating vector..." warning, please just ignore it!
# As long as your RStudio doesn't crash, it is successful.
exprMat <- as.matrix(GetAssayData(cd4_obj, layer = "counts"))

# 2. Extract Meta data
cellInfo <- cd4_obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "seurat_clusters", "detailed_type")] # Remember to add your cell type column

# 3. Extract UMAP coordinates (Optional, but recommended)
umap_coords <- Embeddings(cd4_obj, "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# 4. Generate the Loom file
# This step will generate the file in the current working directory
build_loom(
  file.name = "cd4_scenic_input.loom",
  dgem = exprMat,
  title = "CD4_T_Cells",
  default.embedding = umap_coords,
  default.embedding.name = "UMAP"
)

# 5. Check if the file has been generated
if(file.exists("cd4_scenic_input.loom")) {
  message("✅ Success! The file cd4_scenic_input.loom has been generated, size is approximately: ", 
          round(file.size("cd4_scenic_input.loom")/1024/1024, 2), " MB")
}