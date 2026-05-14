import pandas as pd

import scanpy as sc

import drugreflector as dr



# All major classes available at top level

# dr.DrugReflector, dr.SignatureRefinement, dr.compute_vscores_adata, etc.



# Step 1: Load PBMC 3k dataset with cell type annotations

adata = sc.datasets.pbmc3k()  # Unfiltered dataset with more genes

annots = sc.datasets.pbmc3k_processed().obs  # Cell type annotations



# Merge annotations

adata.obs = pd.merge(adata.obs, annots, how='left', left_index=True, right_index=True)



# Step 2: Compute v-scores between two monocyte populations

vscores = dr.compute_vscores_adata(

    adata, 

    group_col='louvain',

    group1_value='CD14+ Monocytes',    # Classical monocytes

    group2_value='FCGR3A+ Monocytes'   # Non-classical monocytes

)



print(f"V-score comparison: {vscores.name}")

print(f"Computed v-scores for {len(vscores)} genes")

print(f"Top upregulated genes in FCGR3A+ vs CD14+ monocytes:")

print(vscores.nlargest(10))



# Step 3: Initialize DrugReflector with model checkpoints

model_paths = [

    'checkpoints/model_fold_0.pt',

    'checkpoints/model_fold_1.pt', 

    'checkpoints/model_fold_2.pt'

]



model = dr.DrugReflector(checkpoint_paths=model_paths)



# Step 4: Make predictions using v-scores

# DrugReflector will automatically preprocess gene names to HGNC format

predictions = model.predict(vscores, n_top=50)

print(f"Prediction results shape: {predictions.shape}")

print(f"Columns: {predictions.columns.names}")



# Access different metrics

print("\nTop 10 predicted compounds by rank:")

rank_col = ('rank', vscores.name)  # Uses informative name: 'louvain:CD14+ Monocytes->FCGR3A+ Monocytes'

print(predictions[rank_col].nsmallest(10))



print("\nTop 10 compounds by probability:")

prob_col = ('prob', vscores.name)  

print(predictions[prob_col].nlargest(10))



print(f"\nAvailable columns: {list(predictions.columns)}")# 1. Pandas Series (single v-score vector)

# The series name will be used as the transition identifier in outputs

vscore_series = pd.Series([1.2, -0.8, 0.5, ...], index=['GENE1', 'GENE2', 'GENE3', ...], 

                         name='treatment:control->drug')

predictions = model.predict(vscore_series)

# Columns will be: ('rank', 'treatment:control->drug'), ('logit', 'treatment:control->drug'), etc.



# 2. Pandas DataFrame (multiple transitions/signatures)

vscores_df = pd.DataFrame({

    'GENE1': [1.2, 0.8],

    'GENE2': [-0.8, 1.1], 

    'GENE3': [0.5, -0.3]

}, index=['treatment_A', 'treatment_B'])

predictions = model.predict(vscores_df)



# 3. AnnData (v-scores in .X)

vscores_adata = AnnData(

    X=vscores_df.values,

    var=pd.DataFrame(index=vscores_df.columns),

    obs=pd.DataFrame(index=vscores_df.index)

)

predictions = model.predict(vscores_adata)# Example of automatic preprocessing

import pandas as pd

from drugreflector import DrugReflector



# Input with mixed gene name formats

mixed_genes = ['tp53', 'EGFR', 'ENSG00000141510.11', 'CDKN1A_at', 'il6.v2']

vscores = pd.Series([1.2, -0.8, 0.5, 2.1, -1.1], index=mixed_genes)



model = DrugReflector(checkpoint_paths=model_paths)

predictions = model.transform(vscores)



# Output shows preprocessing:

# Preprocessing gene names to HGNC format...

# Preprocessed 4/5 gene names for HGNC compatibility

# Examples of changes:

#   tp53 -> TP53

#   ENSG00000141510.11 -> ENSG00000141510

#   CDKN1A_at -> CDKN1A

#   il6.v2 -> IL6import drugreflector as dr

import pandas as pd



# Starting signature (pandas Series with gene names as index)

starting_signature = pd.Series([1.2, -0.8, 0.5, ...], 

                              index=['GENE1', 'GENE2', 'GENE3', ...])



# Initialize signature refinement (available at top level)

refiner = dr.SignatureRefinement(starting_signature)



# Load experimental data (AnnData with compound treatments)

# adata should have:

# - Gene expression data in .X or layers

# - Compound IDs in .obs (e.g., 'compound_id' column)

# - Sample IDs in .obs (e.g., 'sample_id' column) 

refiner.load_counts_data(

    adata, 

    compound_id_obs_col='compound_id',

    sample_id_obs_cols=['sample_id'],

    layer='raw_counts'  # or None to use .X

)



# Load phenotypic readouts

readouts = pd.Series([0.8, -1.2, 0.3, ...], 

                    index=['compound_A', 'compound_B', 'compound_C', ...])

refiner.load_phenotypic_readouts(readouts)



# Compute learned signatures using correlation analysis

refiner.compute_learned_signatures(corr_method='pearson')



# Generate refined signatures (interpolation between starting and learned)

refiner.compute_refined_signatures(

    learning_rate=0.5,      # 0.5 = equal weight to starting and learned

    scale_learned_sig=True  # Scale learned signature to match starting signature std

)



# Access results

refined_signatures = refiner.refined_signatures  # AnnData object

learned_signatures = refiner.learned_signatures   # AnnData object# For multiple experimental conditions, specify signature_id_obs_cols

refiner.load_counts_data(

    adata,

    compound_id_obs_col='compound_id',

    sample_id_obs_cols=['sample_id'],

    signature_id_obs_cols=['treatment_type', 'timepoint'],  # Creates separate signatures

    layer='raw_counts'

)



# This will create one learned/refined signature for each unique combination

# of values in signature_id_obs_colsfrom drugreflector import compute_vscores_adata, compute_vscore_two_groups



# Compute v-scores between two cell populations

vscores = compute_vscores_adata(

    adata, 

    group_col='cell_type',      # Column identifying groups

    group1_value='control',     # Reference group

    group2_value='treatment',   # Comparison group

    layer=None                  # Use .X, or specify layer name

)



# vscores is a pandas Series with gene names as index and informative name

print(f"V-score comparison: {vscores.name}")  # e.g., "cell_type:control->treatment"

print(f"Top upregulated genes:")

print(vscores.nlargest(10))

print(f"Top downregulated genes:")

print(vscores.nsmallest(10))



# For two arrays directly

group1_values = [1.2, 0.8, 1.5, 0.9]  # Reference/control

group2_values = [2.1, 1.9, 2.3, 2.0]  # Treatment/comparison

vscore = compute_vscore_two_groups(group1_values, group2_values)from drugreflector import load_h5ad_file, pseudobulk_adata



# Load H5AD file with preprocessing

adata = load_h5ad_file('data.h5ad')



# Pseudobulk single-cell data

pseudobulked = pseudobulk_adata(

    adata,

    sample_id_obs_cols=['donor_id', 'condition'],  # Columns defining samples

    method='sum'  # or 'mean'

)from drugreflector import compute_vscores



# Use v-scores in existing workflow

transitions = {

    'group_col': 'cell_type',

    'group1_value': 'control',

    'group2_value': 'treatment'

}



vscores_adata = compute_vscores(adata, transitions=transitions)

# Returns AnnData object with v-scores as .X# Example command line usage (if predict.py exists)

python drugreflector/predict.py input.h5ad \

    --model1 checkpoints/model_fold_0.pt \

    --model2 checkpoints/model_fold_1.pt \

    --model3 checkpoints/model_fold_2.pt \

    --output results.csv