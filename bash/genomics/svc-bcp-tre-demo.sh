#!/bin/bash
# ---------------------------------------------
# Title: Svc_BCP TRE Demo
# Description: From: Source/3. Efforts/_Archives/Svc_BCP TRE Demo.md
# ---------------------------------------------

# PLINK質控流程
plink --vcf simulated_data.vcf \
      --maf 0.01 \           # MAF > 1%
      --geno 0.02 \          # Missing rate < 2%
      --hwe 1e-6 \           # HWE p-value > 1e-6
      --mind 0.02 \          # Sample missing < 2%
      --make-bed \
      --out qc_data