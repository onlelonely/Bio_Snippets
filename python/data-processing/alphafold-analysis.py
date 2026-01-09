# ---------------------------------------------
# Title: AlphaFold
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/AI & Structural Biology/AlphaFold.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 在 Google Colab 中運行
!pip install colabfold[alphafold]
from colabfold import *

# 預測單一蛋白質結構
query_sequence = "MKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL"
results = batch_fold([query_sequence], 
                    job_name="my_protein",
                    model_type="alphafold2_ptm")

# --- Part 2 ---
import numpy as np
from Bio.PDB import PDBParser

def analyze_confidence_scores(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    
    confidence_scores = []
    for residue in structure.get_residues():
        for atom in residue:
            if atom.get_name() == 'CA':  # 只取 alpha carbon
                confidence_scores.append(atom.get_bfactor())
    
    return {
        'mean_confidence': np.mean(confidence_scores),
        'high_confidence_ratio': sum(1 for score in confidence_scores if score > 70) / len(confidence_scores)
    }