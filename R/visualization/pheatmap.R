# ---------------------------------------------
# Title: pheatmap
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/pheatmap.md
# ---------------------------------------------

library(pheatmap)

# 假設 `data_matrix` 是你的數據 (rows: genes/pathways, cols: samples)
# 例如 ssGSEA_results 或 基因表現矩陣

# 1. (可選) 對數據進行縮放 (通常對 row 進行)
#    t() 轉置 -> scale() -> t() 轉置回來
scaled_matrix <- t(scale(t(data_matrix)))

# 2. 準備樣本註釋 (Annotation)
#    建立一個 data frame，其 row names 必須對應到 data_matrix 的 column names
sample_groups <- sub("_.*", "", colnames(data_matrix)) # 從 "Ctrl_1" 提取 "Ctrl"
annotation_col <- data.frame(
  Group = factor(sample_groups),
  row.names = colnames(data_matrix)
)

# 3. (可選) 定義註釋的顏色
ann_colors <- list(
  Group = c(Ctrl = "#1f77b4", Exo = "#ff7f0e", IV = "#2ca02c", Sham = "#d62728")
)

# 4. 繪製熱圖
pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,   # 加上樣本註釋
  annotation_colors = ann_colors,    # 指定註釋顏色
  show_rownames = TRUE,              # 顯示行名 (基因/路徑)
  show_colnames = TRUE,              # 顯示列名 (樣本)
  cluster_cols = TRUE,               # 對列進行群集
  cluster_rows = TRUE,               # 對行進行群集
  fontsize_row = 8,                  # 行字體大小
  main = "My Awesome Heatmap"        # 圖標題
)

# 若要存成檔案
# pdf("my_heatmap.pdf", width=10, height=15)
# pheatmap(...)
# dev.off()