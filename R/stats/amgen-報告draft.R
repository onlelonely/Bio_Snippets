# ---------------------------------------------
# Title: Amgen 報告Draft
# Description: From: Source/3. Efforts/Res_Amgen_MTAP_KRAS/_Archives/Amgen 報告Draft.md
# ---------------------------------------------

# --- Table 1: 病患基本特徵表 ---
# 建立標準的 Table 1 格式，比較 MTAP 狀態間的基本特徵差異

create_table1 <- function(cancer_filter = NULL) {
  # 準備資料
  table1_data <- master_df
  
  if (!is.null(cancer_filter)) {
    table1_data <- table1_data %>% filter(submitted_diagnosis == cancer_filter)
  }
  
  # 合併所有相關資訊
  table1_data <- table1_data %>%
    left_join(death_data %>% select(hash_id, Date_of_death), by = "hash_id") %>%
    left_join(access_critical_biomarker %>% select(hash_id, percent_pd_11), by = "hash_id") %>%
    left_join(biomarkers_tmb_msi %>% select(hash_id, tmb_score, msi_status), by = "hash_id") %>%
    left_join(mtap_isoforms %>% select(hash_id, isoform_simplified), by = "hash_id") %>%
    mutate(
      # MTAP 分組
      mtap_group = case_when(
        mtap_status == "WT" ~ "MTAP Wild-type",
        !is.na(isoform_simplified) ~ "MTAP Altered",
        TRUE ~ "Unknown"
      ),
      
      # 年齡分組
      age_at_diagnosis = 2023 - as.numeric(Year_of_Birth),
      age_group = case_when(
        age_at_diagnosis < 50 ~ "<50",
        age_at_diagnosis >= 50 & age_at_diagnosis < 70 ~ "50-69",
        age_at_diagnosis >= 70 ~ "≥70",
        TRUE ~ "Unknown"
      ),
      
      # 性別
      gender_label = case_when(
        Gender == "1" ~ "Male",
        Gender == "2" ~ "Female",
        TRUE ~ "Unknown"
      ),
      
      # 臨床分期
      stage_group = case_when(
        Clinical_stage %in% c("3A", "3B") ~ "Stage III",
        Clinical_stage %in% c("4", "4A") ~ "Stage IV",
        Clinical_stage == "999" ~ "Unknown",
        TRUE ~ "Other"
      ),
      
      # BMI 計算
      bmi = as.numeric(Weight_kg) / (as.numeric(Height_cm) / 100)^2,
      bmi_group = case_when(
        bmi < 18.5 ~ "Underweight",
        bmi >= 18.5 & bmi < 25 ~ "Normal",
        bmi >= 25 & bmi < 30 ~ "Overweight",
        bmi >= 30 ~ "Obese",
        TRUE ~ "Unknown"
      ),
      
      # PD-L1 分組
      pd_l1_value = as.numeric(percent_pd_11),  # 修正欄位名稱
      pd_l1_group = case_when(
        pd_l1_value >= 50 ~ "High (≥50%)",
        pd_l1_value >= 1 ~ "Low (1-49%)",
        pd_l1_value < 1 ~ "Negative (<1%)",
        TRUE ~ "Unknown"
      ),
      
      # TMB 分組
      tmb_group = case_when(
        tmb_score >= 10 ~ "High (≥10)",
        tmb_score < 10 ~ "Low (<10)",
        TRUE ~ "Unknown"
      ),
      
      # MSI 狀態
      msi_group = case_when(
        msi_status == "MSI-High" ~ "MSI-High",
        msi_status == "MS-Stable" ~ "MSS",
        TRUE ~ "Unknown"
      ),
      
      # 死亡狀態
      death_status = ifelse(!is.na(Date_of_death), "Deceased", "Alive"),
      
      # 癌症類型
      cancer_name = case_when(
        Cancer_type == "01" ~ "Non-squamous NSCLC",
        Cancer_type == "02" ~ "Esophageal cancer",
        Cancer_type == "03" ~ "Gastric cancer",
        Cancer_type == "04" ~ "Gallbladder cancer",
        Cancer_type %in% c("05", "07", "08") ~ "Cholangiocarcinoma",
        Cancer_type == "06" ~ "Pancreatic cancer",
        TRUE ~ "Other"
      )
    ) %>%
    filter(mtap_group != "Unknown")
  
  # 計算統計摘要
  create_summary_stats <- function(data, group_var) {
    data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        percentage = round(n / sum(n) * 100, 1),
        summary = paste0(n, " (", percentage, "%)")
      )
  }
  
  # 建立 Table 1 結構
  table1_results <- list()
  
  # 總體樣本數
  table1_results$total <- table1_data %>%
    group_by(mtap_group) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = mtap_group, values_from = n, values_fill = 0)
  
  # 年齡
  table1_results$age <- table1_data %>%
    group_by(mtap_group) %>%
    summarise(
      mean_age = round(mean(age_at_diagnosis, na.rm = TRUE), 1),
      sd_age = round(sd(age_at_diagnosis, na.rm = TRUE), 1),
      median_age = round(median(age_at_diagnosis, na.rm = TRUE), 1),
      age_summary = paste0(mean_age, " ± ", sd_age),
      .groups = 'drop'
    ) %>%
    select(mtap_group, age_summary) %>%
    pivot_wider(names_from = mtap_group, values_from = age_summary)
  
  # 年齡分組
  table1_results$age_group <- table1_data %>%
    count(mtap_group, age_group) %>%
    group_by(mtap_group) %>%
    mutate(
      percentage = round(n / sum(n) * 100, 1),
      summary = paste0(n, " (", percentage, "%)")
    ) %>%
    select(mtap_group, age_group, summary) %>%
    pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  
  # 性別
  table1_results$gender <- table1_data %>%
    count(mtap_group, gender_label) %>%
    group_by(mtap_group) %>%
    mutate(
      percentage = round(n / sum(n) * 100, 1),
      summary = paste0(n, " (", percentage, "%)")
    ) %>%
    select(mtap_group, gender_label, summary) %>%
    pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  
  # 臨床分期
  table1_results$stage <- table1_data %>%
    count(mtap_group, stage_group) %>%
    group_by(mtap_group) %>%
    mutate(
      percentage = round(n / sum(n) * 100, 1),
      summary = paste0(n, " (", percentage, "%)")
    ) %>%
    select(mtap_group, stage_group, summary) %>%
    pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  
  # BMI
  table1_results$bmi <- table1_data %>%
    group_by(mtap_group) %>%
    summarise(
      mean_bmi = round(mean(bmi, na.rm = TRUE), 1),
      sd_bmi = round(sd(bmi, na.rm = TRUE), 1),
      bmi_summary = paste0(mean_bmi, " ± ", sd_bmi),
      .groups = 'drop'
    ) %>%
    select(mtap_group, bmi_summary) %>%
    pivot_wider(names_from = mtap_group, values_from = bmi_summary)
  
  # PD-L1
  table1_results$pd_l1 <- table1_data %>%
    filter(!is.na(pd_l1_value)) %>%
    count(mtap_group, pd_l1_group) %>%
    group_by(mtap_group) %>%
    mutate(
      percentage = round(n / sum(n) * 100, 1),
      summary = paste0(n, " (", percentage, "%)")
    ) %>%
    select(mtap_group, pd_l1_group, summary) %>%
    pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  
  # TMB
  table1_results$tmb <- table1_data %>%
    group_by(mtap_group) %>%
    summarise(
      median_tmb = round(median(tmb_score, na.rm = TRUE), 1),
      q25_tmb = round(quantile(tmb_score, 0.25, na.rm = TRUE), 1),
      q75_tmb = round(quantile(tmb_score, 0.75, na.rm = TRUE), 1),
      tmb_summary = paste0(median_tmb, " (", q25_tmb, "-", q75_tmb, ")"),
      .groups = 'drop'
    ) %>%
    select(mtap_group, tmb_summary) %>%
    pivot_wider(names_from = mtap_group, values_from = tmb_summary)
  
  # MSI
  table1_results$msi <- table1_data %>%
    count(mtap_group, msi_group) %>%
    group_by(mtap_group) %>%
    mutate(
      percentage = round(n / sum(n) * 100, 1),
      summary = paste0(n, " (", percentage, "%)")
    ) %>%
    select(mtap_group, msi_group, summary) %>%
    pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  
  # 癌症類型 (如果分析所有癌別)
  if (is.null(cancer_filter)) {
    table1_results$cancer_type <- table1_data %>%
      count(mtap_group, cancer_name) %>%
      group_by(mtap_group) %>%
      mutate(
        percentage = round(n / sum(n) * 100, 1),
        summary = paste0(n, " (", percentage, "%)")
      ) %>%
      select(mtap_group, cancer_name, summary) %>%
      pivot_wider(names_from = mtap_group, values_from = summary, values_fill = "0 (0%)")
  }
  
  return(table1_results)
}

# --- 統計檢驗函數 ---
perform_statistical_tests <- function(data, cancer_filter = NULL) {
  if (!is.null(cancer_filter)) {
    data <- data %>% filter(submitted_diagnosis == cancer_filter)
  }
  
  # 準備資料
  test_data <- data %>%
    left_join(death_data %>% select(hash_id, Date_of_death), by = "hash_id") %>%
    left_join(access_critical_biomarker %>% select(hash_id, percent_pd_11), by = "hash_id") %>%
    left_join(biomarkers_tmb_msi %>% select(hash_id, tmb_score, msi_status), by = "hash_id") %>%
    left_join(mtap_isoforms %>% select(hash_id, isoform_simplified), by = "hash_id") %>%
    mutate(
      mtap_group = case_when(
        mtap_status == "WT" ~ "MTAP Wild-type",
        !is.na(isoform_simplified) ~ "MTAP Altered",
        TRUE ~ "Unknown"
      ),
      age_at_diagnosis = 2023 - as.numeric(Year_of_Birth),
      gender_label = case_when(
        Gender == "1" ~ "Male",
        Gender == "2" ~ "Female",
        TRUE ~ "Unknown"
      ),
      pd_l1_value = as.numeric(percent_pd_11),  # 修正欄位名稱
      bmi = as.numeric(Weight_kg) / (as.numeric(Height_cm) / 100)^2
    ) %>%
    filter(mtap_group != "Unknown")
  
  # 進行統計檢驗
  test_results <- list()
  
  # 年齡 (t-test)
  if (length(unique(test_data$mtap_group)) == 2) {
    age_test <- t.test(age_at_diagnosis ~ mtap_group, data = test_data)
    test_results$age_p <- round(age_test$p.value, 3)
    
    # TMB (Wilcoxon test)
    tmb_test <- wilcox.test(tmb_score ~ mtap_group, data = test_data)
    test_results$tmb_p <- round(tmb_test$p.value, 3)
    
    # BMI (t-test)
    bmi_test <- t.test(bmi ~ mtap_group, data = test_data)
    test_results$bmi_p <- round(bmi_test$p.value, 3)
    
    # 性別 (Fisher's exact test)
    gender_table <- table(test_data$mtap_group, test_data$gender_label)
    if (all(dim(gender_table) == c(2, 2))) {
      gender_test <- fisher.test(gender_table)
      test_results$gender_p <- round(gender_test$p.value, 3)
    }
    
    # PD-L1 (Chi-square test)
    pd_l1_data <- test_data %>% filter(!is.na(pd_l1_value))
    if (nrow(pd_l1_data) > 0) {
      pd_l1_table <- table(pd_l1_data$mtap_group, pd_l1_data$pd_l1_value >= 50)
      if (all(dim(pd_l1_table) == c(2, 2))) {
        pd_l1_test <- fisher.test(pd_l1_table)
        test_results$pd_l1_p <- round(pd_l1_test$p.value, 3)
      }
    }
    
    # MSI (Fisher's exact test)
    msi_table <- table(test_data$mtap_group, test_data$msi_status)
    if (all(dim(msi_table) == c(2, 2))) {
      msi_test <- fisher.test(msi_table)
      test_results$msi_p <- round(msi_test$p.value, 3)
    }
  }
  
  return(test_results)
}

# --- 生成各癌別的 Table 1 ---
cancer_types_table1 <- c("Non-squamous NSCLC", "Pancreatic cancer", "Gastric cancer")

for (cancer_type in cancer_types_table1) {
  cat("\n=== Table 1 for", cancer_type, "===\n")
  
  # 生成 Table 1
  table1_result <- create_table1(cancer_type)
  
  # 顯示結果
  cat("Sample Size:\n")
  print(table1_result$total)
  
  cat("\nAge (mean ± SD):\n")
  print(table1_result$age)
  
  cat("\nAge Groups:\n")
  print(table1_result$age_group)
  
  cat("\nGender:\n")
  print(table1_result$gender)
  
  cat("\nClinical Stage:\n")
  print(table1_result$stage)
  
  cat("\nBMI (mean ± SD):\n")
  print(table1_result$bmi)
  
  cat("\nPD-L1 Expression:\n")
  print(table1_result$pd_l1)
  
  cat("\nTMB (median, IQR):\n")
  print(table1_result$tmb)
  
  cat("\nMSI Status:\n")
  print(table1_result$msi)
  
  # 統計檢驗
  test_results <- perform_statistical_tests(master_df, cancer_type)
  cat("\nStatistical Test Results:\n")
  print(test_results)
  
  cat("\n", rep("=", 50), "\n")
}

# --- 生成總體 Table 1 ---
cat("\n=== Overall Table 1 (All Cancer Types) ===\n")
overall_table1 <- create_table1()
overall_tests <- perform_statistical_tests(master_df)

cat("Sample Size:\n")
print(overall_table1$total)

cat("\nCancer Types:\n")
print(overall_table1$cancer_type)

cat("\nAge (mean ± SD):\n")
print(overall_table1$age)

cat("\nGender:\n")
print(overall_table1$gender)

cat("\nTMB (median, IQR):\n")
print(overall_table1$tmb)

cat("\nStatistical Test Results:\n")
print(overall_tests)