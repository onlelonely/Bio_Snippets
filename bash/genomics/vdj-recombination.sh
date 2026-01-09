#!/bin/bash
# ---------------------------------------------
# Title: V(D)J recombination
# Description: From: Source/1. Atlas/🛡️ Immunology & Vaccines/Immune Cells/V(D)J recombination.md
# ---------------------------------------------

# 假設已安裝IGoR
# -m: 指定物種與基因座 (e.g., human TRB)
# -d: 輸入的TCR-seq數據 (通常是已經過初步處理的CDR3序列列表)
# -o: 輸出模型參數檔案

igor -m human TRB -d my_repertoire_data.tsv -o inferred_model.txt