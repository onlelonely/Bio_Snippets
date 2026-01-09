#!/bin/bash
# ---------------------------------------------
# Title: 環境設定指南
# Description: From: Source/3. Efforts/Res_Ab_developability _ML/Archives/環境設定指南.md
# ---------------------------------------------

# Create comprehensive bioinformatics environment
conda create -n h2a-pipe python=3.11 -y
conda activate h2a-pipe

# Core scientific computing
conda install -c conda-forge \
  numpy=1.24.3 \
  pandas=2.0.3 \
  scipy=1.11.1 \
  matplotlib=3.7.2 \
  seaborn=0.12.2 \
  jupyter=1.0.0 \
  ipywidgets=8.0.7 -y

# Bioinformatics specific
conda install -c bioconda \
  biopython=1.81 \
  pysam=0.21.0 \
  bedtools=2.31.0 \
  samtools=1.17 \
  blast=2.14.0 -y

# Machine learning
pip install \
  scikit-learn==1.3.0 \
  xgboost==1.7.6 \
  lightgbm==4.0.0 \
  optuna==3.3.0 \
  shap==0.42.1

# Deep learning
pip install \
  torch==2.0.1 \
  transformers==4.31.0 \
  fair-esm==2.0.0 \
  pytorch-lightning==2.0.6

# Structural biology
conda install -c conda-forge \
  pymol-open-source=2.5.0 \
  mdanalysis=2.5.0 \
  prody=2.4.0 -y

# Additional tools
pip install \
  rdkit==2023.3.2 \
  logomaker==0.8 \
  adjustText==0.8 \
  umap-learn==0.5.3 \
  h5py==3.9.0 \
  tables==3.8.0 \
  dvc==3.15.0 \
  mlflow==2.5.0