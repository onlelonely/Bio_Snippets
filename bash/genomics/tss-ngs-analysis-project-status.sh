#!/bin/bash
# ---------------------------------------------
# Title: TSS NGS Analysis Project Status
# Description: From: Source/3. Efforts/_Archives/Res_TSS NGS Analysis/TSS NGS Analysis Project Status.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 1. 執行過濾後的 BLAST 分析 (推薦)

./run_blast_clean.sh

  

# 2. 如需要下載完整參考序列 (可選)

wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

mkdir -p reference && mv Homo_sapiens.GRCh38.cdna.all.fa reference/

# --- Part 2 ---
# 如果過濾後比對仍不理想，可嘗試:

./run_gtf_blast.sh              # GTF 直接比對

./setup_blast_analysis.sh       # 查看其他工具選項