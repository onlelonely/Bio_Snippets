#!/bin/bash
# ---------------------------------------------
# Title: Human Leukocyte Antigen
# Description: From: Source/1. Atlas/🛡️ Immunology & Vaccines/Immune Cells/Human Leukocyte Antigen.md
# ---------------------------------------------

# 假設已安裝OptiType及所需依賴 (e.g., RazerS3, HDF5)
# -i: 輸入的FASTQ檔案 (R1, R2)
# --dna: 標示為DNA-Seq數據
# -o: 輸出目錄

optitype -i sample_R1.fastq.gz sample_R2.fastq.gz --dna -o /path/to/output/