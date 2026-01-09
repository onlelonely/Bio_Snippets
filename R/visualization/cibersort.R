# ---------------------------------------------
# Title: CIBERSORT
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/CIBERSORT.md
# ---------------------------------------------

# 安裝套件
# remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)
library(readr)

# 載入基因表現數據 (以範例數據為例)
# 假設 `gene_expression_matrix.txt` 是一個基因ID為列名，樣本ID為行名的矩陣
# gene_expression_matrix <- read_tsv("path/to/your/gene_expression_matrix.txt")

# 執行 CIBERSORT
# 注意：要使用 CIBERSORT，你需要從官方網站下載原始碼和 LM22 特徵矩陣
# 並將 `CIBERSORT.R` 和 `LM22.txt` 放在工作目錄下
# set_cibersort_binary("path/to/CIBERSORT.R")
# set_cibersort_mat("path/to/LM22.txt")

# res_cibersort <- deconvolute(gene_expression_matrix, "cibersort")

# 視覺化結果
# library(ggplot2)
# res_cibersort %>% 
#   gather(sample, fraction, -cell_type) %>%
#   ggplot(aes(x=sample, y=fraction, fill=cell_type)) + 
#     geom_bar(stat="identity") + 
#     coord_flip() + 
#     scale_fill_brewer(palette="Paired")