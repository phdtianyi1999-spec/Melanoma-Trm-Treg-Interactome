import scanpy as sc
import cell2location
import numpy as np
import os

# Read the reference file you just generated
adata_ref = sc.read_h5ad("D:/Bioinfo_Temp/reference_melanoma_C2L.h5ad")

# Preparation: ensure gene names are unique, and filter out genes with extremely low expression
adata_ref.var_names_make_unique()

# cell2location recommends retaining only genes expressed in a certain proportion of cells in the reference
from cell2location.utils.filtering import filter_genes

# Remove non_archival parameter
# cell_count_cutoff: genes expressed in at least 5 cells
# cell_percentage_cutoff2: genes expressed in at least 3% of cells
selected_genes = filter_genes(
    adata_ref, 
    cell_count_cutoff=5, 
    cell_percentage_cutoff2=0.03
)

# Confirm the number of genes after filtering
print(f"Number of genes after filtering: {len(selected_genes)}")

# Execute filtering and copy the object
adata_ref = adata_ref[:, selected_genes].copy()



# ================= Configuration Area =================
# Root directory for R export files
base_dir = "D:/Data_us/Exported_Files/"

# Your sample list
sample_list = ["spatial_T13", "spatial_T5", "spatial_T14"]

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scanpy as sc
import pandas as pd
import sys
import warnings
import matplotlib as mpl
import torch
import scvi
import cell2location
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
from cell2location import run_colocation
# ================= 1. Global Configuration Area =================
base_dir = "D:/Data_us/Exported_Files/"
output_dir = "D:/Data_us/spatial_dataus_out\\"  # Location to save the final h5ad

# Define samples and their specific parameters
# The logic here: T13 uses the perfect parameters you determined.
# T5 and T14 temporarily use T13's parameters (if slice positions differ, you only need to modify shift_x/y here)
sample_configs = {
    "spatial_T13": {
        "shift_x": 1800, 
        "shift_y": 2000, 
        "scale_adj": 0.95
    },
    "spatial_T5": {
        "shift_x": 1000,  # If T5 is misaligned, modify this
        "shift_y": 2000, 
        "scale_adj": 0.94 # Scale factor is usually the same for the same batch of data, keep unchanged
    },
    "spatial_T14": {
        "shift_x": 1000,  # If T14 is misaligned, modify this
        "shift_y": 2000, 
        "scale_adj": 0.97
    }
}

# ================= 2. Core Assembly Function (Keep unchanged) =================
def assemble_spatial_final(sample_id, data_dir, shift_x, shift_y, scale_adj):
    print(f"\n>>> Processing sample: {sample_id} ...")
    print(f"    Applying parameters: Shift(X={shift_x}, Y={shift_y}), ScaleAdj={scale_adj}")
    
    # 1. Read data
    mtx_path = os.path.join(data_dir, f"{sample_id}_matrix.mtx")
    if not os.path.exists(mtx_path):
        print(f"!!! Skipping: File not found {mtx_path}")
        return None

    adata = sc.read_mtx(mtx_path).T 
    genes = pd.read_csv(os.path.join(data_dir, f"{sample_id}_features.csv"), header=0, index_col=0)
    barcodes = pd.read_csv(os.path.join(data_dir, f"{sample_id}_barcodes.csv"), header=0, index_col=0)
    adata.obs_names = barcodes.iloc[:, 0].values
    adata.var_names = genes.iloc[:, 0].values
    adata.var_names_make_unique()
    
    # 2. Read and correct coordinates
    coords_path = os.path.join(data_dir, f"{sample_id}_tissue_positions.csv")
    coords_df = pd.read_csv(coords_path, index_col=0)
    common_cells = adata.obs_names.intersection(coords_df.index)
    adata = adata[common_cells].copy()
    coords_df = coords_df.loc[common_cells]
    
    if 'x' in coords_df.columns and 'y' in coords_df.columns:
        spatial_coords = coords_df[['x', 'y']].values
    else:
        spatial_coords = coords_df.iloc[:, 0:2].values

    # Apply Shift
    spatial_coords[:, 0] = spatial_coords[:, 0] - shift_x
    spatial_coords[:, 1] = spatial_coords[:, 1] - shift_y
    adata.obsm['spatial'] = spatial_coords.astype(float)
    
    # 3. Read image
    img_path = os.path.join(data_dir, f"{sample_id}_hires_image.png")
    image = plt.imread(img_path)
    img_h, img_w, _ = image.shape
    
    # 4. Calculate Scale factor
    max_coord_x = np.max(adata.obsm['spatial'][:, 0])
    max_coord_y = np.max(adata.obsm['spatial'][:, 1])
    
    if max_coord_x <= 0 or max_coord_y <= 0:
        print("!!! Warning: Shift value is too large causing coordinates to become negative, please reduce Shift value!")
    
    base_scale_x = img_w / max_coord_x
    base_scale_y = img_h / max_coord_y
    final_scale_factor = min(base_scale_x, base_scale_y) * scale_adj
    
    # 5. Assemble uns
    adata.uns['spatial'] = {
        sample_id: {
            'images': {'hires': image},
            'scalefactors': {
                'tissue_hires_scalef': final_scale_factor,
                'spot_diameter_fullres': 30,
            }
        }
    }
    adata.obs['sample'] = sample_id
    return adata

# ================= 3. Batch Loop Execution =================
processed_adatas = {}

for sample_id, params in sample_configs.items():
    try:
        # Call assembly function
        adata_fixed = assemble_spatial_final(
            sample_id, 
            base_dir, 
            shift_x=params['shift_x'], 
            shift_y=params['shift_y'], 
            scale_adj=params['scale_adj']
        )
        
        if adata_fixed is not None:
            # 1. Save the repaired file
            save_path = os.path.join(output_dir, f"{sample_id}_rebuilt_final.h5ad")
            adata_fixed.write(save_path)
            print(f"    Saved: {save_path}")
            
            # 2. Immediate plot check (generate a plot to help you confirm if T5 and T14 are also aligned)
            plt.figure(figsize=(6, 6))
            sc.pl.spatial(
                adata_fixed, 
                library_id=sample_id,
                color=None,
                spot_size=None,
                title=f"{sample_id} Check (Shift={params['shift_x']}/{params['shift_y']})",
                show=False,
                alpha_img=0.8
            )
            plt.show()
            
            # Store in dictionary for backup
            processed_adatas[sample_id] = adata_fixed
            
    except Exception as e:
        print(f"!!! Error occurred while processing {sample_id}: {e}")

print("\n>>> All samples processed! Please check the generated images above. <<<")
print(">>> If T5 or T14 is misaligned, please modify the shift value in sample_configs and rerun.")


#================================== Start Comments ==================================
import cell2location
import pandas as pd
import os

# ================= 1. Setup Model (Reference) =================
# Ensure labels_key is the cell type column you prepared in R (e.g., 'C2L_final_annotation')
# batch_key: If your scRNA-seq data contains multiple batches/samples, enter the corresponding column name (e.g., "batch");
# If none or already integrated, enter None.
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    labels_key='C2L_final_annotation', 
    batch_key=None 
)

# ================= 2. Initialize and Train Model =================
print(">>> Starting training single-cell reference model...")
mod_ref = cell2location.models.RegressionModel(adata_ref)

# Use GPU accelerated training
# max_epochs=250 is usually enough for the model to converge
# === Optimization 1: Enable Tensor Cores acceleration (for RTX 5080) ===
torch.set_float32_matmul_precision('medium') 

# ================= Train Reference Model =================
print(">>> Starting training single-cell reference model (RTX 5080 Turbo Mode)...")

# Assuming mod_ref is already initialized
# mod_ref = cell2location.models.RegressionModel(adata_ref)
mod_ref.train(max_epochs=2500, batch_size=2500, accelerator="gpu")

# ================= 3. Export Expression Signatures =================
print(">>> Training completed, exporting cell type signatures...")
adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})

# Extract core results: average gene expression profile for each cell type
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg']
else:
    inf_aver = adata_ref.uns['mod']['post_sample_means']['project_col_means']

# This inf_aver is the "dictionary" we will pass to the next step
# Column names are like: "means_per_cluster_mu_fg_[[CD8+ Trm]]"
# We clean up the column names, removing prefixes and keeping only the cell type name
inf_aver.columns = [col.replace("means_per_cluster_mu_fg_[[", "").replace("]]", "") for col in inf_aver.columns]

print(f"Successfully extracted signatures, containing cell types: {inf_aver.columns.tolist()}")

# Perform spatial mapping
import scanpy as sc
import numpy as np

# ================= 1. Configure Paths and Samples =================
# Your specified new path
spatial_dir = "D:/Data_us/spatial_dataus_out/"
# Your three sample names (ensure filename is sample_rebuilt_final.h5ad format, or your previously saved format)
sample_names = ["spatial_T13", "spatial_T5", "spatial_T14"]

# ================= 2. Load and Merge Spatial Data =================
adatas_list = []
print(f">>> Loading spatial data from {spatial_dir}...")

for s in sample_names:
    # Assuming your previously saved filename suffix is _rebuilt_final.h5ad
    # If not, please modify the f-string here
    fname = f"{s}_rebuilt_final.h5ad" 
    fpath = os.path.join(spatial_dir, fname)
    
    if os.path.exists(fpath):
        print(f"    Loading: {fname}")
        ad = sc.read_h5ad(fpath)
        ad.obs['sample'] = s # Ensure sample column exists
        adatas_list.append(ad)
    else:
        print(f"!!! Warning: File not found {fpath}")

# Merge
if len(adatas_list) > 0:
    adata_vis = sc.concat(adatas_list, keys=sample_names, join='outer', index_unique='-')
    
    # *** Key Fix: Restore uns['spatial'] for plotting ***
    adata_vis.uns['spatial'] = {}
    for i, ad in enumerate(adatas_list):
        if 'spatial' in ad.uns:
            adata_vis.uns['spatial'].update(ad.uns['spatial'])
else:
    raise FileNotFoundError("No spatial files loaded, please check path and filenames!")

# ================= 3. Intersect with Reference Signatures =================
# Find shared genes between single-cell signatures and spatial data
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print(f"Number of shared genes: {len(intersect)}")

# Keep only shared genes
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# Ensure Counts are integers
from scipy.sparse import issparse
if issparse(adata_vis.X):
    adata_vis.X.data = np.round(adata_vis.X.data)
else:
    adata_vis.X = np.round(adata_vis.X)

# ================= 4. Train Spatial Mapping Model =================
print(">>> Starting training spatial mapping model (RTX 5080)...")

# Set prior: each Spot has about 30 cells (melanoma)
cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, 
    batch_key="sample" 
)

mod_spatial = cell2location.models.Cell2location(
    adata_vis, 
    cell_state_df=inf_aver, 
    N_cells_per_location=30, 
    detection_alpha=20
)

mod_spatial.train(
    max_epochs=3000, 
    batch_size=None, # Full batch training at once is fastest
    train_size=1,
    accelerator="gpu"
)

# ================= 5. Export Results and Visualization =================
print(">>> Exporting prediction results...")
adata_vis = mod_spatial.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod_spatial.adata.n_obs}
)

# Save final large results
save_path = os.path.join(spatial_dir, "spatial_melanoma_analysis_completed.h5ad")
adata_vis.write(save_path)
print(f"Analysis complete! Results saved to: {save_path}")

# View all predicted cell type names
print(">>> Cell type names included in the results are:")
all_columns = adata_vis.obsm['q05_cell_abundance_w_sf'].columns.tolist()
print(all_columns)

# Plot to show Trm and TH1 Treg
import scanpy as sc
import matplotlib.pyplot as plt
import copy
import os

# ================= 1. Configure Parameters =================
# Must point to the directory where you saved "rebuilt_final" files previously
rebuilt_dir = "D:/Data_us/spatial_dataus_out\\" 
sample_names = ["spatial_T13", "spatial_T5", "spatial_T14"]
target_cells = ['CD8+ Trm', 'Treg Th1 Program']

# Visual fine-tuning parameters
SPOT_SIZE = 8       # Spot size (previously too small, now increased to 3.8)
SPOT_ALPHA = 0.85     # Spot opacity (0-1), higher is more solid
IMG_ALPHA = 0.6       # Background image transparency

# ================= 2. Core Fix: Force Restore Alignment Info =================
print(">>> Forcing restore of alignment parameters from original files...")

# Ensure uns['spatial'] container exists
if 'spatial' not in adata_vis.uns:
    adata_vis.uns['spatial'] = {}

for s in sample_names:
    rebuilt_path = os.path.join(rebuilt_dir, f"{s}_rebuilt_final.h5ad")
    
    if os.path.exists(rebuilt_path):
        try:
            # 1. Read the old "perfectly aligned" file
            # Note: We only read its uns, not the huge X, so it's very fast
            ad_fixed = sc.read_h5ad(rebuilt_path)
            
            # 2. Force overwrite current object's image info
            # This step will copy over the correct image + scale factors
            adata_vis.uns['spatial'][s] = copy.deepcopy(ad_fixed.uns['spatial'][s])
            
            print(f"    ✅ Successfully restored perfect alignment parameters for {s}!")
        except Exception as e:
            print(f"    ❌ Failed to restore {s}: {e}")
    else:
        print(f"    ⚠️ Warning: File not found {rebuilt_path}, cannot fix alignment.")

# ================= 3. Plotting Code A: Single-Cell Distribution =================
print("\n>>> Starting plotting single-cell type distribution...")

for s in sample_names:
    # Extract subset
    subset = adata_vis[adata_vis.obs['sample'] == s]
    
    print(f"Plotting sample: {s}")
    sc.pl.spatial(
        subset, 
        library_id=s, 
        color=target_cells, 
        cmap='magma', 
        
        # --- Corrected parameters ---
        size=SPOT_SIZE,      # Larger spots
        alpha=SPOT_ALPHA,    # Clearer colors
        alpha_img=IMG_ALPHA, # Fainter background
        img_key="hires",     # Force use of high-res image
        # ------------------
        
        vmin=0, 
        vmax='p99', # Automatically enhance contrast
        title=[f"{s}: {c}" for c in target_cells],
        show=True
    )

# ================= 4. Plotting Code B: Co-localization Analysis =================
print("\n>>> Starting plotting Co-localization...")

for s in sample_names:
    # Must copy, otherwise temporary variables will pollute original data
    subset = adata_vis[adata_vis.obs['sample'] == s].copy()
    
    # Ensure subset has the recently fixed image info again
    # (scanpy slicing sometimes loses uns, double insurance here)
    subset.uns['spatial'] = {s: adata_vis.uns['spatial'][s]}
    
    # Calculate co-localization score (product method)
    subset.obs['Co_localization'] = (
        subset.obs['CD8+ Trm'] * subset.obs['Treg Th1 Program']
    )
    
    sc.pl.spatial(
        subset,
        library_id=s,
        color='Co_localization',
        cmap='viridis', # Use Viridis (yellow-green-blue) for higher distinguishability
        
        # --- Corrected parameters ---
        size=SPOT_SIZE,
        alpha=SPOT_ALPHA,
        alpha_img=IMG_ALPHA,
        img_key="hires",
        # ------------------
        
        title=f"{s}: Trm & Treg Co-localization",
        show=True
    )

# Return to R Analysis
# ================= Configure Export Path =================
export_dir = "D:/Data_us/spatial_dataus_out/"
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

print(">>> Preparing to export data to R...")

# 1. Extract predicted cell abundance (q05 is the most robust statistic)
# Copy first to avoid modifying original object
df_export = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()

# 2. Clean column names (remove those annoying prefixes)
# Change 'q05cell_abundance_w_sf_means_per_cluster_mu_fg_CD8+ Trm' to 'CD8+ Trm'
prefix = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
df_export.columns = [col.replace(prefix, "") for col in df_export.columns]

# 3. Add spatial coordinates (crucial for plotting)
# Ensure coordinate names match Seurat habits
df_export['x'] = adata_vis.obsm['spatial'][:, 0]
df_export['y'] = adata_vis.obsm['spatial'][:, 1]

# 4. Add sample info
df_export['sample'] = adata_vis.obs['sample']

# 5. Export as CSV
csv_path = os.path.join(export_dir, "cell2location_results_for_R.csv")
df_export.to_csv(csv_path)

print(f"✅ Export successful! File located at: {csv_path}")
print("   Contains: Cell abundance matrix + coordinates (x,y) + sample origin")