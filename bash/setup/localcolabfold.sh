#!/bin/bash
# ---------------------------------------------
# Title: LocalColabFold
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/AI & Structural Biology/LocalColabFold.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
Git clone https://GitHub.com/YoshitakaMo/localcolabfold.git 
cd localcolabfold
./install_colabbatch_linux.sh

export PATH="/home/cc1296/localcolabfold/localcolabfold/colabfold-conda/bin:$PATH"

# --- Part 2 ---
colabfold_batch --model-type alphafold2_ptm --num-recycle 3 --use-GPU-relax input.fasta results_ptm_recycle3


# --num-recycle 3: 預設的循環次數 
# --model-type AlphaFold2_ptm: 使用 AlphaFold2 模型並產出 pTM score 
# --use-gpu-relax: 使用 GPU 進行結構鬆弛