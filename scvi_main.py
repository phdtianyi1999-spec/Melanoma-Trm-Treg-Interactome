import scanpy as sc
import scvi
import pandas as pd
import scipy.io
import os
import numpy as np

#================= Set Path, don't forget! =========================================
# 1. Check current path (Current Working Directory)
current_path = os.getcwd()
print(f"📍 Your current path is:\n{current_path}")
print("-" * 50)

# 2. Modify current path (Change Directory)
# Change the path below to the folder you really want to go to, such as the directory where you just read the data
target_path = r"D:\Data_us\data\combined_data" 

try:
    os.chdir(target_path)
    print(f"🚀 Success! Current path has been switched to:\n{os.getcwd()}")
    
    # Check again if your results are here
    if os.path.exists("scvi_export_for_r"):
        print("🎉 Found the output folder 'scvi_export_for_r' in the new path!")
    else:
        print("⚠️ Note: The output folder was not found in the new path, it might have been saved elsewhere earlier.")
        
except FileNotFoundError:
    print(f"❌ Error: Cannot find this path {target_path}, please check spelling.")

# Set core number limit
scvi.settings.num_threads = 8 

# ==========================================
# 1. Read Data (Load Data) - [Modified Part]
# ==========================================
# Note: It's best to add an 'r' prefix for Windows paths, or uniformly use forward slashes
data_dir = r"D:\Data_us\data\combined_data" 

print(f"Reading data: {data_dir} ...")

# 1.1 Read sparse matrix
# Scanpy read_mtx usually reads as (genes x cells), needs transpose .T to become (cells x genes)
print("Loading matrix.mtx ...")
adata = sc.read_mtx(os.path.join(data_dir, "matrix.mtx")).T 

# 1.2 Read gene names (features.tsv) and assign to adata.var_names
print("Loading features.tsv ...")
genes = pd.read_csv(os.path.join(data_dir, "features.tsv"), header=None, sep='\t')
# Assuming the first column is gene names
adata.var_names = genes[0].values
# If there are duplicate gene names, this step is crucial, otherwise errors will occur later
adata.var_names_make_unique()

# 1.3 Read cell names (barcodes.tsv) and assign to adata.obs_names
print("Loading barcodes.tsv ...")
barcodes = pd.read_csv(os.path.join(data_dir, "barcodes.tsv"), header=None, sep='\t')
adata.obs_names = barcodes[0].values

# 1.4 Read Metadata
print("Loading metadata.csv ...")
meta_path = os.path.join(data_dir, "metadata.csv")
metadata = pd.read_csv(meta_path, index_col=0)

# Assign Metadata to adata
# Use .loc to ensure the order is consistent with adata.obs_names
adata.obs = metadata.loc[adata.obs_names]

# Check Batch column
BATCH_KEY = 'orig.ident' 

print(f"Data loading complete! Dimensions: {adata.shape}")
print(f"Batch column name used: {BATCH_KEY}")
# Check if this column exists to prevent errors
if BATCH_KEY in adata.obs.columns:
    print("Included batches:", adata.obs[BATCH_KEY].unique())
else:
    print(f"Warning: Cannot find column name '{BATCH_KEY}', please check adata.obs.columns")

# ==========================================
# 2. Preprocessing
# ==========================================
# Filter low-quality cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate mitochondria
# Note: If your gene names are all uppercase (like HUMAN), use startswith('MT-')
# If it's mouse (like Mouse), it's usually 'mt-', please adjust according to the actual situation
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Simple filtering
adata = adata[adata.obs.pct_counts_mt < 20, :]

# --- Key step: Backup Raw Counts ---
adata.layers["counts"] = adata.X.copy()

# Normalize and logarithmize (Only used for finding highly variable genes)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes (HVG)
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=3000, 
    subset=False, 
    layer="counts", 
    flavor="seurat_v3",
    batch_key=BATCH_KEY
)

# ==========================================
# 3. scVI Integration
# ==========================================
# Create HVG subset for training
adata_hvg = adata[:, adata.var['highly_variable']].copy()
# ==========================================
# Optimization step: Convert matrix to CSR format to accelerate training
# ==========================================
# ==========================================
# Then run your original setup_anndata
# ==========================================
from scipy.sparse import csr_matrix

print("Checking data format...")
# Force convert counts layer to CSR (Compressed Sparse Row) format
# This step is the key to making scVI run fast!
if not isinstance(adata.layers["counts"], csr_matrix):
    print("Detected non-CSR format, converting (this may take 1-2 minutes)...")
    adata.layers["counts"] = csr_matrix(adata.layers["counts"])
    print("Conversion complete! It is now in CSR format.")
else:
    print("Data is already in CSR format, no conversion needed.")

scvi.model.SCVI.setup_anndata(
    adata_hvg, 
    layer="counts", 
    batch_key=BATCH_KEY
)

# Initialize model
model = scvi.model.SCVI(adata_hvg, n_latent=10, gene_likelihood="nb")

# Train model
import torch
print(torch.cuda.is_available())
print("Starting to train scVI model...")
model.train(max_epochs=None, accelerator="auto", devices=1)

# Extract features
latent = model.get_latent_representation()
adata.obsm['X_scVI'] = latent

print("scVI integration complete.")

# ==========================================
# 4. Export data for R (Seurat) usage
# ==========================================
output_dir = "scvi_export_for_r"
os.makedirs(output_dir, exist_ok=True)

# Export Embeddings
pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs_names).to_csv(
    os.path.join(output_dir, "scvi_embeddings.csv")
)

# Export Metadata
adata.obs.to_csv(os.path.join(output_dir, "metadata.csv"))

print(f"Results exported to {output_dir}.")