# ---------------------------------------------
# Title: MSigDB
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/MSigDB.md
# ---------------------------------------------

library(msigdbr)
library(dplyr) # 方便進行管線操作

# 1. 查詢可用的物種和類別
# msigdbr_species()
# msigdbr_collections()

# 2. 取得小鼠 (Mus musculus) 的 Hallmark (H) 基因集
H_gene_sets_df <- msigdbr(species = "Mus musculus", category = "H")

# 3. 取得小鼠的 C2 (Curated) 中的 KEGG 子類別基因集
C2_KEGG_df <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")

# 4. 將 data frame 轉換為 GSEA/GSVA 所需的 list 格式
#    列表的名稱 (name) 是基因集名稱 (gs_name)
#    列表的內容 (value) 是該基因集包含的基因符號 (gene_symbol)
gene_sets <- C2_KEGG_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

# `gene_sets` 就可以直接用於 `GSVA::gsva()` 等分析函數