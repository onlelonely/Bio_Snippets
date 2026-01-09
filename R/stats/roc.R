# ---------------------------------------------
# Title: ROC
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/ROC.md
# ---------------------------------------------

# 假設 `predictions` 是模型輸出的機率分數向量
# `labels` 是真實的標籤向量 (0或1)

library(pROC)

# 創建ROC物件
roc_obj <- roc(labels, predictions)

# 計算AUC值
auc_value <- auc(roc_obj)
print(paste("AUC:", round(auc_value, 4)))

# 繪製ROC曲線
plot(roc_obj, main="ROC Curve", print.auc=TRUE)

# 尋找最佳閾值 (例如，使用Youden's J statistic)
best_threshold <- coords(roc_obj, "best", ret="threshold")
print(paste("Best Threshold:", best_threshold))