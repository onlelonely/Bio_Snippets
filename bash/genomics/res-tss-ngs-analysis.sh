#!/bin/bash
# ---------------------------------------------
# Title: Res_TSS NGS Analysis
# Description: From: Source/3. Efforts/_Archives/Res_TSS NGS Analysis/Res_TSS NGS Analysis.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
chmod +x process.sh
nohup ./process.sh > process.log 2>&1 &

# --- Part 2 ---
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gunzip gencode.v47.annotation.gtf.gz

# htseq-count只算UTR
htseq-count -f bam -r pos -s no -t UTR -i gene_name 3F1-3B1.sorted.bam gencode.v47.annotation.gtf > utr_counts_htseq_3F1_3B1.txt


# featureCounts只算UTR
featureCounts -p -a gencode.v47.annotation.gtf -t UTR -g gene_name -o utr_counts_3F1_3B1.txt 3F1-3B1.sorted.bam