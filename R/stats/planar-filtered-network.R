# ---------------------------------------------
# Title: Planar Filtered Network
# Description: From: Source/1. Atlas/🧬 Genomics & Molecular Bio/Sequencing Technologies/Planar Filtered Network.md
# ---------------------------------------------

library(MEGENA)
library(igraph)

# 假設 `datExpr` 是一個基因表現矩陣 (rows=genes, cols=samples)

# 1. 計算基因間的相關性
# MEGENA::calculate.correlation 是一個方便的函數
ijw <- calculate.correlation(datExpr, doPerm = 10)

# 2. 提取邊列表 (edgelist)
edgelist <- ijw[, 1:3]
colnames(edgelist) <- c("row", "col", "weight")

# 3. 計算 PFN
# 這一步會迭代加邊並檢查平面性
pfn_edgelist <- calculate.PFN(edgelist, doPar = TRUE, num.cores = 4)

# 4. 創建 igraph 物件以供後續分析
g <- graph.data.frame(pfn_edgelist, directed = FALSE)

# 現在 g 就是一個 PFN，可以進行社群檢測、視覺化等分析