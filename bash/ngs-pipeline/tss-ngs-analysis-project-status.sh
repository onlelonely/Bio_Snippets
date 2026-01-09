#!/bin/bash
# ---------------------------------------------
# Title: TSS NGS Analysis Project Status
# Description: From: Source/3. Efforts/_Archives/Res_TSS NGS Analysis/TSS NGS Analysis Project Status.md
# ---------------------------------------------

# 3. 完成 Approach 1 過濾 (如需要)

python src/artifact_filtering_optimized.py --input-dir processed_fastq --output-dir Results/filtered_fastq --analysis-file Results/analysis_results/high_frequency_analysis_detailed.csv --results-dir Results/analysis_results

  

# 4. 交叉驗證 (當有兩種方法結果時)

python src/cross_validation_analysis.py --approach1-results Results/analysis_results/high_frequency_analysis_detailed.csv --approach2-results Results/blast_clean_alignment/blast_clean_best_matches.csv --output-dir Results/validation_results

  

# 5. 生成最終數據集

python src/generate_final_dataset.py --input-dir processed_fastq --validation-results Results/validation_results/consensus_classification.csv --output-dir Results/final_cleaned_data

  

# 6. 生成綜合報告

python src/generate_comprehensive_report.py --results-dir Results/analysis_results --validation-dir Results/validation_results --final-dir Results/final_cleaned_data --output-dir Results/final_report