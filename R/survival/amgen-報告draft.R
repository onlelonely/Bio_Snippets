# ---------------------------------------------
# Title: Amgen 報告Draft
# Description: From: Source/3. Efforts/Res_Amgen_MTAP_KRAS/_Archives/Amgen 報告Draft.md
# ---------------------------------------------

# --- 4.1: 資料準備 - 建立存活分析資料集 ---
# 需要載入follow-ups資料
follow_ups <- read_csv("simulated_Clinical_Follow-ups.csv")

# 建立存活分析資料集
prepare_survival_data <- function(cancer_filter = NULL) {
  survival_df <- master_df
  
  if (!is.null(cancer_filter)) {
    survival_df <- survival_df %>% filter(submitted_diagnosis == cancer_filter)
  }
  
  survival_df <- survival_df %>%
    # 合併事件時間與追蹤時間
    left_join(death_data %>% select(hash_id, Date_of_death), by = "hash_id") %>%
    left_join(follow_ups %>% select(hash_id, last_followup_date), by = "hash_id") %>%
    # 合併 MTAP 亞型
    left_join(mtap_isoforms %>% select(hash_id, isoform_simplified), by = "hash_id") %>%
    # 建立 MTAP 分組
    mutate(
      mtap_group = case_when(
        mtap_status == "WT" ~ "MTAP Wild-type",
        !is.na(isoform_simplified) ~ "MTAP Altered",
        TRUE ~ "Unknown"
      )
    ) %>%
    # 計算存活時間和狀態 (假設日期以天數表示)
    mutate(
      # 將字串轉換為數值
      death_days = as.numeric(Date_of_death),
      followup_days = as.numeric(last_followup_date),
      
      # 計算存活時間（以月為單位）
      OS_time = ifelse(!is.na(death_days), death_days, followup_days) / 30.4,
      OS_status = ifelse(!is.na(death_days), 1, 0), # 1=死亡事件, 0=設限
      
      # 計算診斷時年齡
      age_at_diagnosis = 2023 - as.numeric(Year_of_Birth),
      
      # 性別編碼
      gender_factor = factor(Gender, levels = c("1", "2"), labels = c("Male", "Female"))
    ) %>%
    filter(OS_time > 0, !is.na(OS_time)) # 確保存活時間為正且非缺失
  
  return(survival_df)
}

# --- 4.2: 繪製 Kaplan-Meier 曲線 ---
create_km_plot <- function(survival_data, cancer_name = NULL) {
  # 過濾掉未知組別
  plot_data <- survival_data %>%
    filter(mtap_group != "Unknown")
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # 建立 survival object
  surv_obj <- Surv(time = plot_data$OS_time, event = plot_data$OS_status)
  
  # 擬合模型
  fit <- survfit(surv_obj ~ mtap_group, data = plot_data)
  
  # 計算統計檢驗
  surv_diff <- survdiff(surv_obj ~ mtap_group, data = plot_data)
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  # 視覺化
  plot_title <- ifelse(is.null(cancer_name), 
                       "Overall Survival by MTAP Status",
                       paste("Overall Survival by MTAP Status -", cancer_name))
  
  p <- ggsurvplot(
    fit,
    data = plot_data,
    pval = TRUE,
    risk.table = TRUE,
    legend.title = "MTAP Status",
    title = plot_title,
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    palette = c("MTAP Wild-type" = "#66C2A5", "MTAP Altered" = "#FC8D62"),
    risk.table.height = 0.3,
    tables.theme = theme_void(),
    ggtheme = theme_minimal()
  )
  
  return(p)
}

# --- 4.3: 執行 Cox 比例風險模型 ---
perform_cox_analysis <- function(survival_data, cancer_name = NULL) {
  # 過濾掉未知組別
  cox_data <- survival_data %>%
    filter(mtap_group != "Unknown", !is.na(age_at_diagnosis), !is.na(gender_factor))
  
  if (nrow(cox_data) < 10) {
    return(NULL)
  }
  
  # 建立單變量模型
  cox_univariate <- CoxPH(
    Surv(OS_time, OS_status) ~ mtap_group,
    data = cox_data
  )
  
  # 建立多變量模型
  cox_multivariate <- coxph(
    Surv(OS_time, OS_status) ~ mtap_group + age_at_diagnosis + gender_factor,
    data = cox_data
  )
  
  # 提取結果
  univariate_summary <- summary(cox_univariate)
  multivariate_summary <- summary(cox_multivariate)
  
  # 創建結果表
  results_table <- tibble(
    model = c("Univariate", "Multivariate"),
    hr = c(univariate_summary$conf.int[1], multivariate_summary$conf.int[1]),
    ci_lower = c(univariate_summary$conf.int[3], multivariate_summary$conf.int[3]),
    ci_upper = c(univariate_summary$conf.int[4], multivariate_summary$conf.int[4]),
    p_value = c(univariate_summary$coefficients[5], multivariate_summary$coefficients[1, 5])
  )
  
  return(list(
    univariate = cox_univariate,
    multivariate = cox_multivariate,
    results_table = results_table
  ))
}

# --- 4.4: 分析各癌別的存活情況 ---
cancer_types_fig4 <- c("Non-squamous NSCLC", "Pancreatic cancer", "Gastric cancer")
survival_results <- list()

for (cancer_type in cancer_types_fig4) {
  # 準備資料
  survival_data <- prepare_survival_data(cancer_type)
  
  if (nrow(survival_data) > 0) {
    # 繪製 Kaplan-Meier 曲線
    km_plot <- create_km_plot(survival_data, cancer_type)
    if (!is.null(km_plot)) {
      print(km_plot)
    }
    
    # 執行 Cox 分析
    cox_results <- perform_cox_analysis(survival_data, cancer_type)
    if (!is.null(cox_results)) {
      survival_resultscancer_type <- cox_results
      
      # 顯示結果
      cat("\n", cancer_type, "Cox Regression Results:\n")
      print(cox_results$results_table)
    }
  }
}

# --- 4.5: 統計摘要 ---
# 合併所有存活分析結果
if (length(survival_results) > 0) {
  all_cox_results <- map_dfr(names(survival_results), ~{
    results <- survival_results.x$results_table
    results$cancer_type <- .x
    return(results)
  })
  
  print("Summary of Cox Regression Results Across Cancer Types:")
  print(all_cox_results)
}

# --- 4.6: 繪製森林圖 ---
# 如果有足夠的資料，繪製森林圖
if (length(survival_results) >= 2) {
  forest_data <- all_cox_results %>%
    filter(model == "Univariate") %>%
    mutate(
      cancer_type = factor(cancer_type, levels = cancer_types_fig4)
    )
  
  if (nrow(forest_data) > 0) {
    forest_plot <- ggplot(forest_data, aes(x = cancer_type, y = hr)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Hazard Ratios for MTAP Alteration Across Cancer Types",
        x = "Cancer Type",
        y = "Hazard Ratio (95% CI)",
        caption = "Reference: MTAP Wild-type"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_log10()
    
    print(forest_plot)
  }
}

# --- 4.7: 存活時間中位數比較 ---
survival_summary <- map_dfr(cancer_types_fig4, ~{
  survival_data <- prepare_survival_data(.x)
  
  if (nrow(survival_data) > 0) {
    summary_stats <- survival_data %>%
      filter(mtap_group != "Unknown") %>%
      group_by(mtap_group) %>%
      summarise(
        n = n(),
        events = sum(OS_status),
        median_survival = median(OS_time, na.rm = TRUE),
        q25_survival = quantile(OS_time, 0.25, na.rm = TRUE),
        q75_survival = quantile(OS_time, 0.75, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(cancer_type = .x)
    
    return(summary_stats)
  }
})

if (nrow(survival_summary) > 0) {
  print("Survival Time Summary by Cancer Type and MTAP Status:")
  print(survival_summary)
}