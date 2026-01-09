# ---------------------------------------------
# Title: ssGSEA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/ssGSEA.md
# ---------------------------------------------

library(GSVA)
library(msigdbr) # 通常與 msigdbr 搭配使用以獲取基因集

# 假設 `filtered_matrix` 是你的基因表現矩陣 (row: gene symbols, col: samples)
# 假設 `gene_sets` 是你的基因集列表 (list format)

# 1. 執行 ssGSEA 分析
# GSVA::gsva 函數會自動根據參數判斷執行 ssGSEA, GSVA, 或 Z-score
ssGSEA_results <- GSVA::gsva(
  filtered_matrix,
  gene_sets,
  method = "ssgsea", # 明確指定使用 ssGSEA 方法
  verbose = TRUE
)

# 舊版或另一種寫法是使用 ssgseaParam 物件
# ssgsea_params <- GSVA::ssgseaParam(
#   filtered_matrix,
#   gene_sets
# )
# ssGSEA_results <- GSVA::gsva(ssgsea_params)

# 結果 `ssGSEA_results` 是一個矩陣
# rows 是基因集 (pathways)，columns 是樣本
# 值是每個樣本在每個路徑上的 ssGSEA 富集分數