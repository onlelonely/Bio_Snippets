# ---------------------------------------------
# Title: RMA 標準化
# Description: From: Source/1. Atlas/🧬 Genomics & Molecular Bio/Sequencing Technologies/RMA 標準化.md
# ---------------------------------------------

# 首先，需要讀取所有的 .CEL 原始數據檔案
celfiles <- list.files(".", full = TRUE, pattern="*.CEL")
celData <- read.celfiles(celfiles)

# 接著，對讀入的數據執行 RMA
# 這個函數會自動完成上述的背景校正、標準化和匯總步驟
eset <- oligo::rma(celData)

# 結果 eset 會是一個 ExpressionSet 物件，包含了標準化後的基因表現矩陣