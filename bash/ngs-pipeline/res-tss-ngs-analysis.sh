#!/bin/bash
# ---------------------------------------------
# Title: Res_TSS NGS Analysis
# Description: From: Source/3. Efforts/_Archives/Res_TSS NGS Analysis/Res_TSS NGS Analysis.md (4 blocks)
# ---------------------------------------------

# --- Part 1 ---
#!/bin/bash

# Define samples
samples=("3F2-3B2" "3F3-3B3" "3F4-3B4")

# Loop over each sample
for sample in "${samples[@]}"; do
  echo "Processing $sample..."

  # Cutadapt
  cutadapt -u 28 \
    -o ${sample}_R1_trimmed.fastq.gz \
    -p ${sample}_R2_trimmed.fastq.gz \
    NGS01-${sample}_22HNWKLT4_R1.fastq.gz \
    NGS01-${sample}_22HNWKLT4_R2.fastq.gz

  # Alignment and sorting
  bowtie2 -x hg38_index \
    -1 ${sample}_R1_trimmed.fastq.gz \
    -2 ${sample}_R2_trimmed.fastq.gz | \
  samtools view -bS - | \
  samtools sort -o ${sample}.sorted.bam

  # Index
  samtools index ${sample}.sorted.bam

  # HTSeq-count
  htseq-count -f bam -R pos -s no -t UTR -i gene_name \
    ${sample}.sorted.bam gencode.v47.annotation.gtf > utr_counts_htseq_${sample}.txt

  # featureCounts
  featureCounts -p -a gencode.v47.annotation.gtf -t UTR -g gene_name \
    -o utr_counts_${sample}.txt ${sample}.sorted.bam

  echo "$sample done."
done

echo "All samples processed."

# --- Part 2 ---
cutadapt -u 28 -o 3F1-3B1_R1_trimmed.fastq.gz -p 3F1-3B1_R2_trimmed.fastq.gz NGS01-3F1-3B1_22HNWKLT4_R1.fastq.gz NGS01-3F1-3B1_22HNWKLT4_R2.fastq.gz

# --- Part 3 ---
zcat 3F1-3B1_R1_trimmed.fastq.gz | awk 'NR % 4 == 2 {print length($0)}' | sort| uniq -c

# --- Part 4 ---
# Build alignment index
bowtie2-build hg38.fa.gz hg38_index

# Align, sort, convert to bam
bowtie2 -x hg38_index -1 3F1-3B1_R1_trimmed.fastq.gz -2 3F1-3B1_R2_trimmed.fastq.gz | \
samtools view -bS - | \
samtools sort -o 3F1-3B1.sorted.bam && \
samtools index 3F1-3B1.sorted.bam