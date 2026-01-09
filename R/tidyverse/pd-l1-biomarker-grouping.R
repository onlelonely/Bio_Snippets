# ---------------------------------------------
# Title: Bug_Report
# Description: From: Source/3. Efforts/Res_Amgen_MTAP_KRAS/_Archives/Test_results_simulationData/Bug_Report.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 統一使用正確的欄位名稱
biomarkers_pd_l1 <- access_critical_biomarker %>%
  mutate(
    pd_l1_value = as.numeric(percent_pd_11),  # 修正欄位名稱
    # ... 其他代碼
  )

# --- Part 2 ---
# 在Figure 3和Table 1中修正
biomarkers_pd_l1 <- access_critical_biomarker %>%
  mutate(
    pd_l1_value = as.numeric(percent_pd_11),  # 修正欄位名稱
    pd_l1_group = case_when(
      pd_l1_value >= 50 ~ "High (≥50%)",
      pd_l1_value >= 1  ~ "Low (1-49%)",
      pd_l1_value < 1   ~ "Negative (<1%)",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(hash_id, pd_l1_group, pd_l1_value)