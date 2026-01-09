# ---------------------------------------------
# Title: Hierarchical Clustering
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/Hierarchical Clustering.md
# ---------------------------------------------

# 假設 `data_matrix` 是一個基因表現矩陣 (rows=genes, cols=samples)

library(pheatmap)

# 使用 pheatmap 自動進行階層式分群並繪製熱圖
# pheatmap 內部會調用 hclust
pheatmap(
  data_matrix,
  scale = "row",  # 對基因進行Z-score標準化，以觀察相對表現
  clustering_distance_rows = "correlation", # 對基因使用相關性距離
  clustering_distance_cols = "euclidean",   # 對樣本使用歐氏距離
  clustering_method = "ward.D2",          # 使用Ward法進行連結
  main = "Hierarchical Clustering of Gene Expression"
)