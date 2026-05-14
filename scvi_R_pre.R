## Change system error messages to English
Sys.setenv(LANGUAGE = "en")
## Disable conversion to factors
options(stringsAsFactors = FALSE)
## Clear environment
rm(list=ls())
library(Seurat)
library(Matrix)
library(data.table)
getwd()
setwd("D:/Data_us")
# ==========================================
# 1. Organize and merge data
# ==========================================
# Assuming scRNAlist is a list containing 19 Seurat objects
# We first merge the data in the list
# Note: add.cell.ids is used to prevent duplicate cell names
sce.list.merged <- merge(x = scRNAlist[[1]], 
                         y = scRNAlist[2:length(scRNAlist)], 
                         add.cell.ids = names(scRNAlist), 
                         project = "List_19_Samples")

# Then merge this combined object with your sce.all_raw (8 cases)
# Also add prefixes to prevent barcode conflicts
sce.combined <- merge(x = sce.all_raw, 
                      y = sce.list.merged, 
                      add.cell.ids = c("Batch8", "Batch19"), # You can change the prefixes here yourself
                      project = "Combined_All")
rm(sce.all_raw,sce.list.merged,scRNAlist)
# ==========================================
# 2. Prepare Metadata (Very important)
# ==========================================
# scVI needs to know what the 'batch' is for each row.
# Check if sce.combined$orig.ident can represent your 27 samples?
# If yes, use orig.ident; if not, you need to manually correct it here.
head(sce.combined$orig.ident) 

# ==========================================
# 3. Export to a format readable by Python
# ==========================================
out_dir <- "./data/combined_data"
dir.create(out_dir, recursive = TRUE)

# 3.1 Export sparse matrix (Counts)
# Note: Be sure to export counts (Raw Data), not data (Normalized)
counts_mat <- GetAssayData(sce.combined, slot = "counts")
writeMM(counts_mat, file = file.path(out_dir, "matrix.mtx"))

# 3.2 Export gene names (Features)
write.table(data.frame(rownames(counts_mat), rownames(counts_mat)), 
            file = file.path(out_dir, "features.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3.3 Export cell names (Barcodes)
write.table(colnames(counts_mat), 
            file = file.path(out_dir, "barcodes.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3.4 Export Metadata (includes batch info)
write.csv(sce.combined@meta.data, file = file.path(out_dir, "metadata.csv"))

print("R data export complete! You can now go to Python to run scVI.")
# Save the combined data
save_file <- "sce_combined_backup.rds"
saveRDS(sce.combined, file = save_file, compress = FALSE)