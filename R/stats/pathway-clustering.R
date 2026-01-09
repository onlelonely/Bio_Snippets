# ---------------------------------------------
# Title: Pathway Clustering
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/Pathway Clustering.md
# ---------------------------------------------

# 假設 `ssGSEA_results` 是 ssGSEA 的輸出矩陣 (rows: pathways, cols: samples)

# 1. 計算路徑間的距離並進行層次聚類
#    dist() 計算距離矩陣，hclust() 進行聚類
hc <- hclust(dist(ssGSEA_results))

# 2. (可選) 視覺化樹狀圖，以幫助決定切割點
# pdf("pathway_dendrogram.pdf")
# plot(hc, labels=FALSE, main="Pathway Clustering Dendrogram")
# abline(h=0.5, col="red", lty=2) # 畫一條參考線
# dev.off()

# 3. 切割樹狀圖以獲得群集分配
#    可以基於高度 (h) 或指定的群集數目 (k)
#    `cluster_height` 的值需要根據樹狀圖的視覺結果來調整
cluster_height <- 0.15
cutree_result <- cutree(hc, h = cluster_height)
# table(cutree_result) # 查看每個群集的大小

# 4. (應用) 找出與你感興趣的路徑位於同一個群集中的其他路徑
# 假設 `candidate_pathways` 是一個包含你感興趣的路徑名稱的向量
candidate_indices <- which(rownames(ssGSEA_results) %in% candidate_pathways)

# 找到這些候選路徑所在的群集編號
candidate_clusters <- unique(cutree_result[candidate_indices])

# 篩選出所有位於這些目標群集中的路徑
related_indices <- which(cutree_result %in% candidate_clusters)
filtered_ssGSEA_results <- ssGSEA_results[related_indices, ]

# `filtered_ssGSEA_results` 現在包含了與候選路徑功能相關的路徑子集