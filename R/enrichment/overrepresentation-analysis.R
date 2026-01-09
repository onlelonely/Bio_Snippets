# ---------------------------------------------
# Title: Overrepresentation analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/Overrepresentation analysis.md
# ---------------------------------------------

# 安裝並載入套件
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db") # 以人類為例

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. 準備基因列表
# `gene_list` 是一個包含差異表現基因 Entrez ID 的向量
# `background_genes` 是一個包含所有檢測到的基因 Entrez ID 的向量
# gene_list <- c("100", "1000", "1001", ...)
# background_genes <- c("1", "2", "3", ..., "100", "101", ...)

# 2. 執行 GO ORA
go_ora <- enrichGO(gene         = gene_list,
                  universe     = background_genes,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = 'ENTREZID',
                  ont          = "BP", # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# 3. 視覺化結果
barplot(go_ora, showCategory=15)
dotplot(go_ora, showCategory=15)

# 4. 執行 KEGG ORA
kegg_ora <- enrichKEGG(gene         = gene_list,
                     universe     = background_genes,
                     organism     = 'hsa', # Homo sapiens
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

# 視覺化 KEGG 結果
barplot(kegg_ora, showCategory=15)