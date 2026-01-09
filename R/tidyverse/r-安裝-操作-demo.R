# ---------------------------------------------
# Title: R 安裝 & 操作 Demo
# Description: From: Source/3. Efforts/_Archives/Svc_食藥大數據競賽_2025/R 安裝 & 操作 Demo.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 建立迴歸分析用的資料集，選取特定欄位
regression_label <- label %>%
  filter(`連鎖與否` == 1) %>%
  select(
    "現場調製飲料標示",
    "食品業者登錄",
    "食品業者保存來源文件",
    "GHP",
    "分店數",
    "人口密度人平方公里",
    "業者類別分組",
    "人口密度分組2",
    "分店數分組2"
  ) 

# 將變數轉換為因子
regression_label$`現場調製飲料標示` <- as.factor(regression_label$`現場調製飲料標示`)
regression_label$`食品業者登錄` <- factor(regression_label$`食品業者登錄`, ordered = FALSE)
regression_label$`食品業者保存來源文件` <- factor(regression_label$`食品業者保存來源文件`, ordered = FALSE)
regression_label$`GHP` <- factor(regression_label$`GHP`, ordered = FALSE)
regression_label$`業者類別分組` <- factor(regression_label$`業者類別分組`, ordered = FALSE)
regression_label$`人口密度分組2` <- factor(regression_label$`人口密度分組2`, ordered = FALSE)
regression_label$`分店數分組2` <- factor(regression_label$`分店數分組2`, ordered = FALSE)

# 重新設定因子的參考組別
regression_label$`現場調製飲料標示` <- relevel(regression_label$`現場調製飲料標示`, ref = "合格")
regression_label$`業者類別分組` <- relevel(regression_label$`業者類別分組`, ref = "1")
regression_label$`GHP` <- relevel(regression_label$`GHP`, ref = "1")
regression_label$`食品業者登錄` <- relevel(regression_label$`食品業者登錄`, ref = "1")
regression_label$`分店數分組2` <- relevel(regression_label$`分店數分組2`, ref = "2")
regression_label$`人口密度分組2` <- relevel(regression_label$`人口密度分組2`, ref = "2")

# --- Part 2 ---
# 建立檢驗資料的迴歸分析資料集
regression_test <- test %>%
  select(
    "檢驗結果",
    "連鎖與否",
    "食品業者登錄",
    "食品業者保存來源文件",
    "GHP",
    "飲料類別",
    "業者類別分組",
    "人口密度分組",
    "平均氣溫分組",
    "平均濕度分組"
  ) 

# 轉換為因子
regression_test <- regression_test %>%
  mutate(across(everything(), as.factor))

regression_test$`飲料類別` <- relevel(regression_test$`飲料類別`, ref = '配料')