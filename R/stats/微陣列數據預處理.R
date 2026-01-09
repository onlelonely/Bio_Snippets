# ---------------------------------------------
# Title: 微陣列數據預處理
# Description: From: Source/1. Atlas/🧬 Genomics & Molecular Bio/Sequencing Technologies/微陣列數據預處理.md
# ---------------------------------------------

library(affycoretools)
library(pd.clariom.s.mouse) # 範例晶片的註釋包

# 假設 `eset` 是經 oligo::RMA() 處理後的 ExpressionSet 物件

# 1. 註釋 ExpressionSet
#    fData(eset) 將會多出 SYMBOL, GENENAME 等欄位
all.eset <- annotateEset(eset, annotation(eset))

# 2. 探針篩選 (僅保留 TC 探針)
tc_rows <- grep("^TC", rownames(all.eset), value = TRUE)
tc_eset <- all.eset[tc_rows, ]

# 提取表現矩陣和註釋資訊
raw_expr <- exprs(tc_eset)
probe_reference <- as.data.frame(fData(all.eset))

# 3. 探針 ID 轉為基因符號
#    建立一個新的矩陣，並將 rownames 從 Probe ID 換成 Gene Symbol
#    注意：此時可能產生NA或重複的 rownames
eset_final <- raw_expr
indices <- match(rownames(eset_final), probe_reference$PROBEID)
rownames(eset_final) <- probe_reference$SYMBOL[indices]

# 清理缺失或不必要的行/列
eset_final <- eset_final[!is.na(rownames(eset_final)), ] # 移除沒有對應到基因符號的探針

# 4. 處理重複基因 (以保留平均表現量最高者為例)
# 計算每個基因 (row) 的平均表現量
means <- rowMeans(eset_final)

# 依基因名和平均表現量排序 (高至低)
# 這樣能確保在去重複時，會保留下表現量最高的那筆
ordered_matrix <- eset_final[order(rownames(eset_final), -means), ]

# 去除重複的基因名，由於已排序，保留下來的是第一個，也就是表現量最高的
final_matrix <- ordered_matrix[!duplicated(rownames(ordered_matrix)), ]

# `final_matrix` 就是可以用於下游分析的乾淨基因表現矩陣