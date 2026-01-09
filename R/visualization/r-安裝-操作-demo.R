# ---------------------------------------------
# Title: R 安裝 & 操作 Demo
# Description: From: Source/3. Efforts/_Archives/Svc_食藥大數據競賽_2025/R 安裝 & 操作 Demo.md (7 blocks)
# ---------------------------------------------

# --- Part 1 ---
library(readxl)
library(dplyr)
library(car)
library(glmnet)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(tidyr)
library(scales)
library(viridis)
library(patchwork)
library(rms)
library(pROC)
library(RANN)
library(caret)

# --- Part 2 ---
# 建立更美觀的 LASSO 圖表
plot_data <- data.frame(
  lambda = cv_lasso$lambda,
  cvm = cv_lasso$cvm,
  cvup = cv_lasso$cvup,
  cvlo = cv_lasso$cvlo
)

p1 <- ggplot(plot_data, aes(x = log(lambda))) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup), alpha = 0.2, fill = "#2C3E50") +
  geom_line(aes(y = cvm), color = "#2C3E50", size = 1) +
  geom_vline(xintercept = log(cv_lasso$lambda.min), 
             linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = log(cv_lasso$lambda.1se), 
             linetype = "dashed", color = "#27AE60") +
  theme_minimal() +
  labs(title = "LASSO 交叉驗證結果",
       x = "log(Lambda)",
       y = "Mean-Squared Error") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 建立係數路徑圖
coef_path <- as.matrix(cv_lasso$glmnet.fit$beta)
lambda_path <- cv_lasso$glmnet.fit$lambda
coef_data <- data.frame(lambda = rep(lambda_path, each = nrow(coef_path)),
                       coefficient = as.vector(coef_path),
                       variable = rep(rownames(coef_path), times = length(lambda_path)))

p2 <- ggplot(coef_data, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(size = 1) +
  geom_vline(xintercept = log(cv_lasso$lambda.min), 
             linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = log(cv_lasso$lambda.1se), 
             linetype = "dashed", color = "#27AE60") +
  theme_minimal() +
  scale_color_viridis_d(option = "plasma") +
  labs(title = "LASSO 係數路徑圖",
       x = "log(Lambda)", 
       y = "Coefficients",
       color = "變項") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# 組合圖表
p1 + p2 + plot_layout(ncol = 1, heights = c(1, 1.5))

# --- Part 3 ---
# 最終邏輯斯迴歸模型
logistic_model <- glm(
  `現場調製飲料標示` ~ 
    `食品業者登錄` + 
    `GHP` + 
    `人口密度分組2` + 
    `分店數分組2`,
  family = binomial(link = "logit"),
  data = regression_label
)

# 建立漂亮的摘要表格
summary_stats <- summary(logistic_model)
coef_df <- as.data.frame(summary_stats$coefficients)
coef_df$Variable <- rownames(coef_df)
coef_df$OR <- exp(coef_df$Estimate)
ci <- exp(confint(logistic_model))
coef_df$CI_Lower <- ci[,1]
coef_df$CI_Upper <- ci[,2]

# 重新排列欄位
coef_df <- coef_df[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)", "OR", "CI_Lower", "CI_Upper")]

print(coef_df)

# 建立森林圖 (Forest plot)
forest_data <- coef_df[-1,] # 移除截距項
forest_data$Variable <- factor(forest_data$Variable, levels = forest_data$Variable)

ggplot(forest_data, aes(y = Variable, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#2C3E50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "邏輯斯迴歸分析之勝算比及95%信賴區間",
       x = "勝算比 (log scale)",
       y = "變項")

# --- Part 4 ---
# 建立更美觀的 LASSO 圖表
plot_data <- data.frame(
  lambda = cv_lasso_test$lambda,
  cvm = cv_lasso_test$cvm,
  cvup = cv_lasso_test$cvup,
  cvlo = cv_lasso_test$cvlo
)

p1 <- ggplot(plot_data, aes(x = log(lambda))) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup), alpha = 0.2, fill = "#2C3E50") +
  geom_line(aes(y = cvm), color = "#2C3E50", size = 1) +
  geom_vline(xintercept = log(cv_lasso_test$lambda.min), 
             linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = log(cv_lasso_test$lambda.1se), 
             linetype = "dashed", color = "#27AE60") +
  theme_minimal() +
  labs(title = "LASSO 交叉驗證結果",
       x = "log(Lambda)",
       y = "Mean-Squared Error") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 建立係數路徑圖
coef_path <- as.matrix(cv_lasso_test$glmnet.fit$beta)
lambda_path <- cv_lasso_test$glmnet.fit$lambda
coef_data <- data.frame(lambda = rep(lambda_path, each = nrow(coef_path)),
                       coefficient = as.vector(coef_path),
                       variable = rep(rownames(coef_path), times = length(lambda_path)))

p2 <- ggplot(coef_data, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(size = 1) +
  geom_vline(xintercept = log(cv_lasso_test$lambda.min), 
             linetype = "dashed", color = "#E74C3C") +
  geom_vline(xintercept = log(cv_lasso_test$lambda.1se), 
             linetype = "dashed", color = "#27AE60") +
  theme_minimal() +
  scale_color_viridis_d(option = "plasma") +
  labs(title = "LASSO 係數路徑圖",
       x = "log(Lambda)", 
       y = "Coefficients",
       color = "變項") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# 組合圖表
p1 + p2 + plot_layout(ncol = 1, heights = c(1, 1.5))

# --- Part 5 ---
# 最終邏輯斯迴歸模型
logistic_model_test <- glm(
  `檢驗結果` ~ 
    `食品業者登錄` + 
    `GHP` + 
    `飲料類別` + 
    `業者類別分組` + 
    `人口密度分組` + 
    `平均濕度分組`,
  family = binomial(link = "logit"),
  data = regression_test
)

# 建立漂亮的摘要表格
summary_stats <- summary(logistic_model_test)
coef_df <- as.data.frame(summary_stats$coefficients)
coef_df$Variable <- rownames(coef_df)
coef_df$OR <- exp(coef_df$Estimate)
ci <- exp(confint(logistic_model_test))
coef_df$CI_Lower <- ci[,1]
coef_df$CI_Upper <- ci[,2]

# 重新排列欄位
coef_df <- coef_df[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)", "OR", "CI_Lower", "CI_Upper")]

print(coef_df)

# 建立森林圖
forest_data <- coef_df[-1,] # 移除截距項
forest_data$Variable <- factor(forest_data$Variable, levels = forest_data$Variable)

ggplot(forest_data, aes(y = Variable, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#2C3E50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "邏輯斯迴歸分析之勝算比及95%信賴區間",
       x = "勝算比 (log scale)",
       y = "變項")

# --- Part 6 ---
# 計算 ROC 曲線的預測機率
predicted_probs <- predict(logistic_model_test, type = "response")
roc_obj <- pROC::roc(regression_test$檢驗結果, predicted_probs)
auc_value <- pROC::auc(roc_obj)

# 建立 ROC 曲線資料框
roc_data <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

# 建立增強的 ROC 圖
ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#2C3E50", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.75, y = 0.25,
           label = sprintf("AUC = %.3f\n95%% CI: %.3f - %.3f",
                         auc_value,
                         pROC::ci(roc_obj)[1],
                         pROC::ci(roc_obj)[3]),
           size = 4, color = "#2C3E50") +
  theme_minimal() +
  labs(title = "ROC曲線",
       x = "偽陽性率 (1 - Specificity)",
       y = "真陽性率 (Sensitivity)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_equal()

# --- Part 7 ---
# 計算預測機率
predicted_probs <- predict(f, type="fitted")

# 建立 ROC 分析用的資料框
roc_data <- data.frame(
  actual = as.numeric(regression_test_eng$test_result) - 1,
  predicted = predicted_probs
)

# 計算 ROC 曲線
roc_obj <- pROC::roc(roc_data$actual, roc_data$predicted)

# 使用 Youden's J 統計量找到最佳閾值
coords <- coords(roc_obj, "best", best.method="youden")
optimal_cutoff <- coords$threshold

# 計算指標的函數
calculate_metrics <- function(actual, predicted, threshold) {
  predicted_class <- ifelse(predicted >= threshold, 1, 0)
  conf_matrix <- table(Actual = actual, Predicted = predicted_class)
  
  # 計算指標
  sensitivity <- conf_matrix[2,2] / sum(conf_matrix[2,])
  specificity <- conf_matrix[1,1] / sum(conf_matrix[1,])
  ppv <- conf_matrix[2,2] / sum(conf_matrix[,2])
  npv <- conf_matrix[1,1] / sum(conf_matrix[,1])
  f1 <- 2 * (ppv * sensitivity) / (ppv + sensitivity)
  
  return(list(
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    NPV = npv,
    F1 = f1,
    Confusion_Matrix = conf_matrix
  ))
}

# 在最佳閾值下計算指標
metrics <- calculate_metrics(roc_data$actual, roc_data$predicted, optimal_cutoff)

# 建立結果表格
metrics_df <- data.frame(
  Metric = c("最佳閾值", "敏感度", "特異度", "陽性預測值", "陰性預測值", "F1分數"),
  Value = c(optimal_cutoff, 
           metrics$Sensitivity,
           metrics$Specificity,
           metrics$PPV,
           metrics$NPV,
           metrics$F1)
)

print(metrics_df)

# 顯示混淆矩陣
print(metrics$Confusion_Matrix)

# 建立不同閾值點的圖表
cutoffs <- seq(0, 1, by = 0.01)
metrics_by_cutoff <- sapply(cutoffs, function(cut) {
  tryCatch({
    m <- calculate_metrics(roc_data$actual, roc_data$predicted, cut)
    c(m$Sensitivity, m$Specificity)
  }, error = function(e) {
    # 如果發生錯誤則回傳 NA 值
    c(NA, NA)
  })
})

cutoff_data <- data.frame(
  Cutoff = cutoffs,
  Sensitivity = metrics_by_cutoff[1,],
  Specificity = metrics_by_cutoff[2,]
)

cutoff_data_long <- tidyr::pivot_longer(cutoff_data, 
                                      cols = c(Sensitivity, Specificity),
                                      names_to = "Metric",
                                      values_to = "Value")

ggplot(cutoff_data_long, aes(x = Cutoff, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = optimal_cutoff, linetype = "dashed", color = "gray50") +
  annotate("text", x = optimal_cutoff, y = 0.25,
           label = sprintf("最佳閾值 = %.3f", optimal_cutoff),
           angle = 90, vjust = -0.5) +
  theme_minimal() +
  labs(title = "不同閾值下的敏感度和特異度",
       x = "閾值",
       y = "數值",
       color = "指標") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_viridis_d(labels = c("敏感度", "特異度")) +
  xlim(0, 0.5)