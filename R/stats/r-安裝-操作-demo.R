# ---------------------------------------------
# Title: R 安裝 & 操作 Demo
# Description: From: Source/3. Efforts/_Archives/Svc_食藥大數據競賽_2025/R 安裝 & 操作 Demo.md (14 blocks)
# ---------------------------------------------

# --- Part 1 ---
install.packages(c(
  "Tidyverse",    # 資料前處理與視覺化
  "data.table",   # 高效能資料處理
  "rstatix",      # 友善的統計檢定
  "glmnet",       # Lasso / Elastic Net
  "survival",     # 生存分析 (CoxPH)
  "survminer",    # 生存曲線視覺化
  "rms",          # Nomogram 與迴歸建模
  "readxl",       # Excel 檔案讀取
  "car",          # 回歸診斷
  "kableExtra",   # 美化表格
  "viridis",      # 色彩配置
  "patchwork",    # 圖形組合
  "pROC",         # ROC 曲線分析
  "RANN",         # K-近鄰
  "caret"         # 機器學習工具
))

# --- Part 2 ---
# 取得當前目錄中所有 Excel 檔案
excel_files <- list.files(path = "./", pattern = "\\.xlsx$", full.names = TRUE)

# 建立空的清單來儲存資料框
data_list <- list()

# 迴圈讀取每個 Excel 檔案
for (file in excel_files) {
    # 提取檔案名稱（不含副檔名）用於命名資料框
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # 讀取 Excel 檔案，跳過前 4 行，使用第 5 行作為欄位名稱
    data <- read_excel(file, skip = 4, col_names = TRUE)
    
    # 將資料框儲存到清單中，以檔案名為鍵值
    data_listfile_name <- data
    
    cat("Loaded:", file_name, "\n")
}

label <- data_list1
test <- data_list2

# --- Part 3 ---
# 建立清單來儲存所有結果
chi_square_results <- list()

# 要分析的預測變數清單
predictor_vars <- c("食品業者登錄", "食品業者保存來源文件", "GHP", "業者類別分組", "人口密度分組2", "分店數分組2")

# 對每個預測變數進行迴圈
for(var in predictor_vars) {
  # 建立列聯表
  cont_table <- table(regression_labelvar, regression_label$`現場調製飲料標示`)
  
  # 執行卡方檢定並處理錯誤
  chi_test <- tryCatch({
    chisq.test(cont_table)
  }, error = function(e) {
    # 如果檢定失敗則回傳 NA 值
    list(p.value = NA, statistic = NA, parameter = NA)
  })
  
  # 儲存結果
  chi_square_resultsvar <- list(
    variable = var,
    p_value = as.numeric(chi_test$p.value),
    chi_squared = as.numeric(chi_test$statistic),
    df = as.numeric(chi_test$parameter),
    significant = !is.na(chi_test$p.value) && chi_test$p.value < 0.05
  )
}

# 將結果轉換為資料框
chi_results_df <- data.frame(
  Variable = sapply(chi_square_results, function(x) x$variable),
  Chi_Squared = sapply(chi_square_results, function(x) ifelse(is.na(x$chi_squared), NA, round(x$chi_squared, 2))),
  DF = sapply(chi_square_results, function(x) x$df),
  P_Value = sapply(chi_square_results, function(x) ifelse(is.na(x$p_value), NA, round(x$p_value, 4))),
  Significant = sapply(chi_square_results, function(x) x$significant)
)

# 顯示結果
print(chi_results_df)

# --- Part 4 ---
# 執行邏輯斯迴歸
logistic_model <- glm(
  `現場調製飲料標示` ~ 
    `食品業者登錄` + 
    `食品業者保存來源文件` + 
    `GHP` + 
    `業者類別分組` + 
    `人口密度分組2` + 
    `分店數分組2`,
  family = binomial(link = "logit"),
  data = regression_label
)

# 顯示摘要統計和勝算比
summary_result <- summary(logistic_model)
print(summary_result)

odds_ratios <- exp(coef(logistic_model))
conf_int <- exp(confint(logistic_model))

or_table <- cbind(
  OR = odds_ratios,
  CI_lower = conf_int[, 1],
  CI_upper = conf_int[, 2],
  p_value = summary_result$coefficients[, 4]
)

print(or_table)

# 計算變異數膨脹因子 (VIF)
vif_values <- car::vif(logistic_model)
print(vif_values)

# --- Part 5 ---
# 移除「食品業者保存來源文件」變數並建立新模型
regression_label <- regression_label %>%
  select(-`食品業者保存來源文件`)

logistic_model <- glm(
  `現場調製飲料標示` ~ 
    `食品業者登錄` + 
    `GHP` + 
    `業者類別分組` + 
    `人口密度分組2` + 
    `分店數分組2`,
  family = binomial(link = "logit"),
  data = regression_label
)

new_result <- summary(logistic_model)
print(new_result)

odds_ratios <- exp(coef(logistic_model))
conf_int <- exp(confint(logistic_model))

or_table <- cbind(
  OR = odds_ratios,
  CI_lower = conf_int[, 1],
  CI_upper = conf_int[, 2],
  p_value = new_result$coefficients[, 4]
)

print(or_table)

# --- Part 6 ---
# 檢查並處理遺漏值
missing_values <- colSums(is.na(regression_label))
print(missing_values)

regression_label <- regression_label[complete.cases(regression_label), ]

# 檢查類別分佈
table_result <- table(regression_label$`現場調製飲料標示`)
print(table_result)
cat("Percentage of class 1:", round(table_result[2]/sum(table_result)*100, 2), "%\n")

# 準備 LASSO 迴歸的資料
y <- as.numeric(regression_label$`現場調製飲料標示`) - 1

x <- model.matrix(~ `食品業者登錄` + 
                   `GHP` + 
                   `業者類別分組` + 
                   `人口密度分組2` + 
                   `分店數分組2` - 1, 
                 data = regression_label)

# 執行交叉驗證
cv_lasso <- cv.glmnet(x, y, 
                      family = "binomial", 
                      alpha = 1,
                      nfolds = 5, 
                      nlambda = 1000)

cat("Optimal lambda (minimum MSE):", cv_lasso$lambda.min, "\n")
cat("Optimal lambda (1SE rule):", cv_lasso$lambda.1se, "\n\n")

# --- Part 7 ---
# 變數重要性的 Bootstrap 分析
n_iterations <- 1000
n_samples <- nrow(x)
coef_matrix <- matrix(0, nrow = n_iterations, ncol = ncol(x))
colnames(coef_matrix) <- colnames(x)

for (i in 1:n_iterations) {
  boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
  x_boot <- x[boot_idx, ]
  y_boot <- y[boot_idx]
  
  cv_fit <- cv.glmnet(x_boot, y_boot, family = "binomial", alpha = 1, nfolds = 5)
  coefs <- coef(cv_fit, s = "lambda.min")[-1]
  coef_matrix[i, ] <- as.numeric(coefs != 0)
}

nonzero_percentage <- colMeans(coef_matrix) * 100
importance_df <- data.frame(
  Variable = colnames(x),
  NonZero_Percentage = nonzero_percentage
)
importance_df <- importance_df[order(importance_df$NonZero_Percentage, decreasing = TRUE), ]
print(importance_df)

# --- Part 8 ---
# 建立清單來儲存所有結果
chi_square_results <- list()

# 要分析的預測變數清單
predictor_vars_test <- c("連鎖與否", "食品業者登錄", "食品業者保存來源文件", "GHP", 
                        "飲料類別", "業者類別分組", "人口密度分組", "平均氣溫分組", "平均濕度分組")

# 對每個預測變數進行迴圈並執行檢驗結果的卡方檢定
for (var in predictor_vars_test) {
  # 建立列聯表
  cont_table <- table(regression_testvar, regression_test$`檢驗結果`)
  
  # 執行卡方檢定
  chi_test <- chisq.test(cont_table)
  
  # 儲存結果
  chi_square_resultsvar <- list(
    variable = var,
    chi_square = chi_test$statistic,
    df = chi_test$parameter,
    p_value = chi_test$p.value,
    significant = chi_test$p.value < 0.05
  )
}

# 將結果轉換為資料框以供顯示
chi_square_df <- data.frame(
  Variable = sapply(chi_square_results, function(x) x$variable),
  Chi_Square = sapply(chi_square_results, function(x) x$chi_square),
  DF = sapply(chi_square_results, function(x) x$df),
  P_Value = sapply(chi_square_results, function(x) x$p_value),
  Significant = sapply(chi_square_results, function(x) x$significant)
)

print(chi_square_df)

# --- Part 9 ---
# 對檢驗資料進行邏輯斯迴歸建模
regression_test$`平均濕度分組` <- relevel(regression_test$`平均濕度分組`, ref = 4)
logistic_model_test <- glm(
  `檢驗結果` ~ 
    `連鎖與否` + 
    `食品業者登錄` + 
    `食品業者保存來源文件` + 
    `GHP` + 
    `飲料類別` + 
    `業者類別分組` + 
    `人口密度分組` + 
    `平均氣溫分組` + 
    `平均濕度分組`,
  family = binomial(link = "logit"),
  data = regression_test
)

# 顯示摘要
summary_test_result <- summary(logistic_model_test)
print(summary_test_result)

# 計算勝算比和信賴區間
test_odds_ratios <- exp(coef(logistic_model_test))
test_conf_int <- exp(confint(logistic_model_test))

# 組合勝算比和信賴區間
test_or_table <- cbind(
  OR = test_odds_ratios,
  CI_lower = test_conf_int[, 1],
  CI_upper = test_conf_int[, 2],
  p_value = summary_test_result$coefficients[, 4]
)

print(format(test_or_table, scientific = FALSE, digits = 5))

# --- Part 10 ---
# VIF 分析
vif_test <- vif(logistic_model_test)
print(vif_test)

# --- Part 11 ---
# 檢查遺漏值
cat("\n===== 檢驗資料遺漏值檢查 =====\n\n")
missing_values_test <- colSums(is.na(regression_test))
print(missing_values_test)

# 移除有遺漏值的列
regression_test <- regression_test[complete.cases(regression_test), ]
cat("\n移除不完整案例後，檢驗資料集有", nrow(regression_test), "個觀測值\n")

# 檢查類別不平衡
table_result_test <- table(regression_test$`檢驗結果`)
print(table_result_test)
cat("檢驗資料中 class 1 的百分比:", round(table_result_test[2]/sum(table_result_test)*100, 2), "%\n")

# --- Part 12 ---
y <- as.numeric(regression_test$`檢驗結果`) - 1

# 建立模型矩陣
x <- model.matrix(~ `連鎖與否` + 
                   `食品業者登錄` + 
                   `食品業者保存來源文件` + 
                   `GHP` + 
                   `飲料類別` + 
                   `業者類別分組` + 
                   `人口密度分組` + 
                   `平均氣溫分組` + 
                   `平均濕度分組` - 1, 
                 data = regression_test)

# 執行交叉驗證
cv_lasso_test <- cv.glmnet(x, y,
  family = "binomial",
  alpha = 1,
  nfolds = 5,
  nlambda = 1000
)

cat("\nOptimal lambda (minimum MSE):", cv_lasso_test$lambda.min, "\n")
cat("Optimal lambda (1SE rule):", cv_lasso_test$lambda.1se, "\n\n")

# --- Part 13 ---
# Bootstrap LASSO
n_iterations <- 1000
n_samples <- nrow(x)
coef_matrix <- matrix(0, nrow = n_iterations, ncol = ncol(x))
colnames(coef_matrix) <- colnames(x)

for (i in 1:n_iterations) {
  boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
  x_boot <- x[boot_idx, ]
  y_boot <- y[boot_idx]
  
  cv_fit <- cv.glmnet(x_boot, y_boot, family = "binomial", alpha = 1, nfolds = 5)
  coefs <- coef(cv_fit, s = "lambda.min")[-1]
  coef_matrix[i, ] <- as.numeric(coefs != 0)
}

# 計算重要性
nonzero_percentage <- colMeans(coef_matrix) * 100
importance_df <- data.frame(
  Variable = colnames(x),
  NonZero_Percentage = nonzero_percentage
)
importance_df <- importance_df[order(importance_df$NonZero_Percentage, decreasing = TRUE), ]
print(importance_df)

# --- Part 14 ---
# 建立英文欄位名稱版本的資料
regression_test_eng <- regression_test
names(regression_test_eng) <- c(
  "test_result",
  "is_chain",
  "business_registration",
  "source_documents",
  "GHP",
  "beverage_type",
  "business_category",
  "population_density",
  "avg_temperature",
  "avg_humidity"
)

# 使用 rms 套件建立列線圖
dd <- datadist(regression_test_eng)
options(datadist="dd")
f <- lrm(test_result ~ 
           business_registration + 
           GHP + 
           beverage_type + 
           business_category + 
           population_density + 
           avg_humidity,
         data=regression_test_eng)
nomogram_p <- nomogram(f, fun=function(x)1/(1+exp(-x)),
                      funlabel="Probability")

# 建立英文到中文名稱的對應
name_mapping <- c(
  "business_registration" = "食品業者登錄",
  "GHP" = "GHP",
  "beverage_type" = "飲料類別", 
  "business_category" = "業者類別分組",
  "population_density" = "人口密度分組",
  "avg_humidity" = "平均濕度分組",
  "total.points" = "total.points",
  "lp" = "lp",
  "Probability" = "Probability"
)

# 在列線圖中以中文名稱替換英文名稱
names(nomogram_p)[match(names(name_mapping), names(nomogram_p))] <- name_mapping[match(names(nomogram_p), names(name_mapping))]

# 建立格式更佳的列線圖
plot(nomogram_p,
     xfrac = 0.2,
     yfrac = 0.2, 
     col.grid = "#2C3E50",
     lwd = 2,
     cex.axis = 0.8,
     cex.var = 0.8,
     cex.lab = 0.8,
     cex.cat = 0.8)