# ---------------------------------------------
# Title: GO
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/GO.md
# ---------------------------------------------

# 安裝套件
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db") # 以人類為例

library(clusterProfiler)
library(org.Hs.eg.db)

# 假設 `gene_list` 是一個包含差異表現基因 Entrez ID 的向量
# gene_list <- c("100", "1000", "1001", ...)

# GO 富集分析 (Over-Representation Analysis)
ego <- enrichGO(gene         = gene_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP", # 可選 "BP", "MF", "CC" 或 "ALL"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

# 視覺化結果
barplot(ego, showCategory=20)
dotplot(ego, showCategory=20)

# 將結果轉換為數據框
ego_df <- as.data.frame(ego)