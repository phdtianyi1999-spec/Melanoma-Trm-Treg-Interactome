# Custom Scripts for: Trm-instructed Treg Reprogramming in Melanoma Immunotherapy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the custom bioinformatics pipeline and analytical scripts used in our study investigating the tumor microenvironment (TME) dynamics during neoadjuvant anti-PD-1 therapy in melanoma. 

Our study integrates single-nucleus RNA sequencing (snRNA-seq), spatial transcriptomics, and an AI-guided drug discovery framework to reveal how CD8+ tissue-resident memory T (Trm) cells orchestrate the phenotypic reprogramming of regulatory T (Treg) cells into a Th1-like state. Furthermore, we validated the prognostic and predictive value of the Trm-Treg interactome across independent external clinical cohorts, translating these findings into a novel combinatorial therapeutic strategy (NSC693868).

The provided scripts cover a cross-platform (R and Python) workflow, including data integration (scVI), high-resolution sub-clustering, trajectory inference, gene regulatory network analysis (pySCENIC), cell-cell communication (NicheNet), spatial deconvolution (cell2location), deep learning-assisted drug screening (DrugReflector), and clinical survival analysis.

## System Requirements

### Hardware Requirements
The scripts require a standard workstation or high-performance computing (HPC) environment. 
* **RAM:** 64+ GB recommended (due to large Seurat/AnnData objects and spatial datasets).
* **CPU:** 8+ cores recommended for parallel processing (especially for pySCENIC GRN inference).
* **GPU:** Strongly recommended (NVIDIA RTX series or better) for Python-based deep learning steps (`scVI`, `cell2location`, and `DrugReflector`).

### Software Requirements
* **OS:** Linux (Ubuntu 20.04+ recommended for server-side pySCENIC) or Windows 10/11.
* **R version:** >= 4.1.0 
* **Python version:** >= 3.8

### Major Dependencies 
**R Packages:**
* `Seurat` (v4.3.0 or v5.0.0), `monocle3`, `nichenetr`, `scDblFinder`, `SeuratDisk`, `msigdbr`.
* Visualization & Stats: `ggplot2`, `patchwork`, `ggpubr`, `pheatmap`, `pROC`, `survival`, `survminer`.

**Python Packages:**
* `scanpy`, `scvi-tools`, `cell2location`, `pyscenic`, `drugreflector`, `torch`, `pandas`, `numpy`.

## Installation Guide

### R Environment Setup
```R
# Install CRAN packages
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "msigdbr", "pals", "ggrepel", "pROC", "survival", "survminer", "pheatmap", "ggpubr"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDblFinder", "AUCell", "ComplexHeatmap"))

# Install GitHub specific packages
devtools::install_github("cole-trapnell-lab/monocle3")
devtools::install_github("saeyslab/nichenetr")
remotes::install_github("mojaveazure/seurat-disk")
Python Environment Setup (via Conda)
Bash
conda create -n melanoma_env python=3.10.19
conda activate melanoma_env
pip install scanpy scvi-tools cell2location pyscenic torch pandas numpy
# For DrugReflector installation, please refer to its official documentation.
Typical installation time: 30 - 60 minutes, depending on network speed and pre-installed dependencies.

Repository Structure & Usage Instructions
The pipeline is modularized into sequential R and Python scripts. We encourage users to follow this cross-platform workflow to reproduce the quantitative results and manuscript figures.

1. Preprocessing and Deep-Learning Data Integration (scVI)
scvi_R_pre.R: Prepares and exports the raw count matrices and metadata from R to Python.

scvi_main.py: Executes deep generative modeling via scvi-tools (SCVI) to correct batch effects across 19 longitudinal samples. Returns a batch-corrected latent representation.

scvi_R_post.R: Imports the scVI latent embeddings back into R, performs doublet removal using scDblFinder, and executes standard Seurat clustering and preliminary cell type annotation.

batch_effect.R: Generates UMAP visualizations to validate the effectiveness of batch-effect correction.

2. High-Resolution T Cell Sub-clustering
Tcells_anno.R: Isolates the pan-T cell compartment, performs sub-clustering, and meticulously annotates distinct subsets (e.g., CD8+ Trm, Th1-like Tregs, CD4+ T_ISG) based on canonical markers.

3. Trajectory Inference
Tcells_Trajectory_Inference.R: Utilizes monocle3 to construct pseudotime trajectories, revealing the phenotypic transition from activated/stem-like Tregs to Th1-like transitional Tregs.

4. Gene Regulatory Network Analysis (SCENIC)
scenic_cd4_r2pyscenic.R: Converts Seurat objects to .loom format for Python.

pyScenic_Server.py: A robust, multiprocessing Python script tailored for Linux servers. Executes the pySCENIC pipeline (GRNBoost2, CTX, AUCell) to infer transcription factor regulatory networks.

scenic_cd4_pyscenic2r.R: Imports the Regulon Activity Scores (AUC) back into R. Calculates Regulon Specificity Scores (RSS) to identify key transcription factors driving Treg fragility.

5. Cell-Cell Communication (NicheNet)
NicheNet_pepiline2.R: Infers intercellular communication networks. Identifies the Trm-derived TGF-β and IFN-γ signaling axes driving Treg phenotypic reprogramming.

6. Spatial Transcriptomics Deconvolution (cell2location)
cell2location_R2python.R: Prepares matched spatial transcriptomics (Visium) and snRNA-seq reference data into .h5ad format.

cell2location_main.py: Executes the Bayesian probabilistic model cell2location in Python to map single-cell signatures onto spatial spots. Validates the in situ spatial co-localization of Trm cells and Th1-like Tregs.

cell2location_py2R.R: Imports the cell2location abundance matrices back into R Seurat spatial objects. Computes spatial co-localization scores, generates custom high-contrast spatial feature plots (e.g., RGB overlay projections), and performs spatial correlation analysis to validate the in situ interaction between CD8+ Trm cells and Th1-like Tregs.

7. AI-Guided Drug Discovery & Pathway Scoring
pathway_score.R: Employs AddModuleScore coupled with msigdbr to quantify specific functional states (IFN-γ response, TGF-β signaling, Treg fragility).

Drugreflector.py: Utilizes the deep learning framework DrugReflector to perform in silico virtual screening, predicting pharmacological compounds (e.g., NSC693868) capable of inducing the Treg phenotypic shift.

Figure4D-4E,6A.R: Extracts global pathway activity scores to generate heatmaps (Fig 4D-4E) and visualizes the top 100 pharmacological candidates predicted by DrugReflector (Fig 6A).

8. Clinical Cohort Validation & Survival Analysis
Figure_7_Riaz.R: Calculates the composite Trm and Th1-like Treg gene signature within external immunotherapy cohorts. Generates violin plots and ROC curves to evaluate predictive efficacy for ICB responsiveness.

Figure7_multi.R: Integrates multiple external clinical datasets, performs survival analysis (Overall Survival and Progression-Free Survival) using the composite signature, and generates Kaplan-Meier curves and multi-ROC comparisons across different baseline models.

Data Availability
Due to patient privacy regulations, our raw sequencing data generated in this study have been deposited in the Genome Sequence Archive (GSA) at the National Genomics Data Center (NGDC), China. These datasets are already publicly accessible to all readers globally under the following accession codes: scRNA-seq data [HRA005837] and spatial transcriptomics data [HRA004456]. 

License
This project is licensed under the MIT License - see the LICENSE file for details.