#!/bin/bash
# ---------------------------------------------
# Title: T-cell
# Description: From: Source/1. Atlas/🛡️ Immunology & Vaccines/Immune Cells/T-cell.md
# ---------------------------------------------

# 步驟1: 使用MiXCR對原始FASTQ檔案進行比對
mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true --report alignment_report.log R1.fastq.gz R2.fastq.gz alignments.vdjca

# 步驟2: 組裝克隆 (Assemble Clones)
mixcr assemble --report assemble_report.log alignments.vdjca clones.clns

# 步驟3: 導出克隆型數據
mixcr exportClones --chains TRB clones.clns clones.tsv