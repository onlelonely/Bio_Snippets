# ---------------------------------------------
# Title: Amgen 報告Draft
# Description: From: Source/3. Efforts/Res_Amgen_MTAP_KRAS/_Archives/Amgen 報告Draft.md (3 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 篩選出 MTAP 相關變異
mtap_variants <- variants_combined %>%
  filter(gene == "MTAP")

# 計算每位病患的變異類型數量
mtap_class_summary <- mtap_variants %>%
  group_by(hash_id) %>%
  summarise(
    alteration_types = list(unique(alteration_type)),
    n_types = n_distinct(alteration_type),
    .groups = 'drop' 
  ) %>%
  mutate(
    mtap_class = case_when(
      n_types > 1 ~ "Multiple",
      map_lgl(alteration_types, ~ "SV" %in% .x) ~ "SV",
      map_lgl(alteration_types, ~ "CNA" %in% .x) ~ "CNA",
      map_lgl(alteration_types, ~ "RE" %in% .x) ~ "RE",
      TRUE ~ "Other"
    )
  )

# 將變異分類資訊與主資料表合併
analysis_fig1 <- master_df %>%
  left_join(mtap_class_summary, by = "hash_id")

# --- Figure 1a ---
# 數一下Top 15 total mutations count
# 先找出所有 MTAP 變異的病患
mtap_patient_ids <- variants_combined %>%
  filter(gene == "MTAP") %>%
  distinct(hash_id) %>%
  pull(hash_id)

# 從所有變異中，篩選出上述病患，並排除 MTAP 基因本身
co_occurring_variants <- variants_combined %>%
  filter(hash_id %in% mtap_patient_ids) %>%
  filter(gene != "MTAP")

# 計算每個基因被多少個 '不同' 病患所擁有
gene_counts <- co_occurring_variants %>%
  group_by(gene) %>%
  summarise(n_patients = n_distinct(hash_id)) %>%
  arrange(desc(n_patients))

# 取出前 15 名
top_15_genes <- gene_counts %>%
  top_n(15, n_patients)

# 計算總 MTAP 病患數，用於計算百分比
total_mtap_patients <- length(mtap_patient_ids)

# 繪圖
ggplot(top_15_genes, aes(x = reorder(gene, n_patients), y = (n_patients / total_mtap_patients) * 100)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 Co-occurring Alterations in MTAP-Altered Tumors",
    x = "Gene",
    y = "Frequency (%)"
  ) +
  theme_classic() +
  # 在 bar 旁邊加上百分比文字
  geom_text(
    aes(label = sprintf("%.1f%%", (n_patients / total_mtap_patients) * 100)),
    hjust = -0.2,
    size = 3.5
  ) +
  scale_y_continuous(limits = c(0, max((top_15_genes$n_patients / total_mtap_patients) * 100) * 1.15))

# --- Figure 1b ---
prevalence_by_cancer <- analysis_fig1 %>%
  group_by(Cancer_type) %>%
  summarise(
    total_cases = n(),
    mtap_altered_cases = sum(!is.na(mtap_class)),
    prevalence_pct = (mtap_altered_cases / total_cases) * 100
  ) %>%
  arrange(desc(prevalence_pct))

ggplot(prevalence_by_cancer, aes(x = reorder(Cancer_type, prevalence_pct), y = prevalence_pct)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Prevalence of MTAP Alterations by Cancer Type",
       x = "Cancer Type", y = "Prevalence (%)") +
  theme_classic()

# --- Figure 1c/d ---
# 註解：原始使用KRAS突變的正規表達式，實際上需要根據真實MTAP變異資料調整
# 暫時改為MTAP變異的次數分配表

# # 從 short variants 表中提取特定 MTAP 突變
# mtap_isoforms <- variant_short %>%
#   filter(gene == "MTAP") %>%
#   select(hash_id, protein_effect) %>%
#   mutate(isoform = str_extract(protein_effect, "G12[A-Z]|G13[A-Z]|Q61[A-Z]")) %>%
#   mutate(isoform_simplified = case_when(
#     isoform %in% c("G12D", "G12V", "G12C", "G13D", "G12R") ~ isoform,
#     !is.na(isoform) ~ "Other",
#     TRUE ~ NA_character_
#   )) %>%
#   filter(!is.na(isoform_simplified))

# MTAP變異次數分配表
mtap_variant_distribution <- variant_short %>%
  filter(gene == "MTAP") %>%
  left_join(master_df %>% select(hash_id, Cancer_type), by = "hash_id") %>%
  # 根據Cancer Type編碼篩選目標癌別並重新分類
  filter(Cancer_type %in% c("01", "02", "03", "04", "05", "06", "07", "08")) %>%
  mutate(cancer_category = case_when(
    Cancer_type == "01" ~ "Non-squamous NSCLC",
    Cancer_type == "02" ~ "Esophageal cancer", 
    Cancer_type == "03" ~ "Gastric cancer",
    Cancer_type == "04" ~ "Gallbladder cancer",
    Cancer_type %in% c("05", "07", "08") ~ "Cholangiocarcinoma",
    Cancer_type == "06" ~ "Pancreatic cancer",
    TRUE ~ "Other"
  )) %>%
  count(cancer_category, protein_effect, sort = TRUE) %>%
  group_by(cancer_category) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# 顯示MTAP變異分佈表
print("MTAP Protein Effect Distribution by Cancer Type:")
mtap_variant_distribution %>%
  arrange(cancer_category, desc(n)) %>%
  select(cancer_category, protein_effect, n, percentage) %>%
  knitr::kable(digits = 1, col.names = c("Cancer Type", "Protein Effect", "Count", "Percentage (%)"))

# 簡化的視覺化 - 顯示前10個最常見的變異
top_mtap_variants <- mtap_variant_distribution %>%
  group_by(protein_effect) %>%
  summarise(total_n = sum(n)) %>%
  top_n(10, total_n) %>%
  pull(protein_effect)

mtap_variant_distribution %>%
  filter(protein_effect %in% top_mtap_variants) %>%
  ggplot(aes(x = cancer_category, y = n, fill = protein_effect)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "MTAP Protein Effect Distribution in Major Cancers",
       x = "Cancer Type", y = "Number of Cases", fill = "MTAP Protein Effect") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

# 暫時建立 mtap_isoforms 變數以供後續分析使用
# 注意：實際分析需要根據真實MTAP變異特性重新定義
mtap_isoforms <- variant_short %>%
  filter(gene == "MTAP") %>%
  select(hash_id, protein_effect) %>%
  mutate(isoform_simplified = "MTAP_altered") # 暫時統一標記為MTAP變異

# --- Part 2 ---
# --- 2.1: 資料準備 - 建立基因突變矩陣 ---
# 分析所有癌別的MTAP共現情況
create_gene_matrix <- function(cancer_filter = NULL) {
  # 篩選特定癌別（如果指定）
  if (!is.null(cancer_filter)) {
    target_patients <- master_df %>%
      filter(submitted_diagnosis == cancer_filter) %>%
      pull(hash_id)
  } else {
    target_patients <- master_df$hash_id
  }
  
  # 建立基因突變矩陣
  gene_matrix <- variants_combined %>%
    filter(hash_id %in% target_patients) %>%
    select(hash_id, gene) %>%
    filter(!is.na(gene)) %>%
    distinct() %>%
    mutate(mutated = 1) %>%
    pivot_wider(names_from = gene, values_from = mutated, values_fill = 0)
  
  # 確保MTAP欄位存在
  if (!"MTAP" %in% colnames(gene_matrix)) {
    gene_matrix$MTAP <- 0
  }
  
  # 更新MTAP狀態
  gene_matrix <- gene_matrix %>%
    mutate(MTAP = ifelse(hash_id %in% mtap_mutated_patients$hash_id, 1, 0)) %>%
    select(hash_id, MTAP, everything())
  
  return(gene_matrix)
}

# --- 2.2: 執行 Fisher's Exact Test ---
perform_cooccurrence_analysis <- function(gene_matrix, min_mutations = 5) {
  co_occurrence_results <- list()
  other_genes <- setdiff(colnames(gene_matrix), c("hash_id", "MTAP"))
  
  for (gene_name in other_genes) {
    # 檢查基因突變頻率
    gene_mut_count <- sum(gene_matrixgene_name)
    if (gene_mut_count < min_mutations) next
    
    # 建立 2x2 contingency table
    contingency_table <- table(
      MTAP_status = gene_matrix$MTAP,
      Other_gene_status = gene_matrixgene_name
    )
    
    # 執行 Fisher's Exact Test
    if(nrow(contingency_table) == 2 && ncol(contingency_table) == 2) {
      test_result <- fisher.test(contingency_table)
      
      # 計算co-occurrence metrics
      both_mutated <- sum(gene_matrix$MTAP == 1 & gene_matrixgene_name == 1)
      mtap_only <- sum(gene_matrix$MTAP == 1 & gene_matrixgene_name == 0)
      gene_only <- sum(gene_matrix$MTAP == 0 & gene_matrixgene_name == 1)
      neither <- sum(gene_matrix$MTAP == 0 & gene_matrixgene_name == 0)
      
      co_occurrence_resultsgene_name <- tibble(
        gene = gene_name,
        p_value = test_result$p.value,
        odds_ratio = test_result$estimate,
        ci_lower = test_result$conf.int[1],
        ci_upper = test_result$conf.int[2],
        both_mutated = both_mutated,
        mtap_only = mtap_only,
        gene_only = gene_only,
        neither = neither,
        gene_mut_freq = gene_mut_count / nrow(gene_matrix)
      )
    }
  }
  
  return(bind_rows(co_occurrence_results))
}

# --- 2.3: 分析各癌別的共現情況 ---
cancer_types <- c("Non-squamous NSCLC", "Pancreatic cancer", "Gastric cancer")
all_results <- list()

for (cancer_type in cancer_types) {
  gene_matrix <- create_gene_matrix(cancer_type)
  results <- perform_cooccurrence_analysis(gene_matrix)
  
  if (nrow(results) > 0) {
    results$cancer_type <- cancer_type
    all_resultscancer_type <- results
  }
}

# 合併結果
combined_results <- bind_rows(all_results)

# --- 2.4: 繪製改良的 Volcano Plot ---
if (nrow(combined_results) > 0) {
  plot_data <- combined_results %>%
    mutate(
      log2_OR = log2(pmax(odds_ratio, 0.01)),  # 避免log(0)
      neg_log10_p = -log10(pmax(p_value, 1e-10)),  # 避免-log10(0)
      significance = case_when(
        p_value < 0.001 ~ "P < 0.001",
        p_value < 0.01 ~ "P < 0.01",
        p_value < 0.05 ~ "P < 0.05",
        TRUE ~ "Not Significant"
      ),
      relationship = case_when(
        log2_OR > 0 & p_value < 0.05 ~ "Co-occurrence",
        log2_OR < 0 & p_value < 0.05 ~ "Mutual Exclusivity",
        TRUE ~ "No Association"
      )
    )
  
  # 按癌別分別繪製
  for (cancer_type in cancer_types) {
    cancer_data <- plot_data %>% filter(cancer_type == !!cancer_type)
    
    if (nrow(cancer_data) > 0) {
      p <- ggplot(cancer_data, aes(x = log2_OR, y = neg_log10_p, color = relationship)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        labs(
          title = paste("MTAP Co-occurrence Analysis -", cancer_type),
          x = "Log2(Odds Ratio)",
          y = "-log10(p-value)",
          color = "Relationship"
        ) +
        theme_minimal() +
        scale_color_manual(values = c(
          "Co-occurrence" = "#E31A1C",
          "Mutual Exclusivity" = "#1F78B4", 
          "No Association" = "#999999"
        ))
      
      # 標示顯著的基因
      significant_genes <- cancer_data %>%
        filter(p_value < 0.05) %>%
        arrange(p_value) %>%
        head(10)
      
      if (nrow(significant_genes) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = significant_genes,
          aes(label = gene),
          size = 3,
          max.overlaps = 10
        )
      }
      
      print(p)
    }
  }
}

# --- 2.5: 統計摘要 ---
if (nrow(combined_results) > 0) {
  summary_stats <- combined_results %>%
    group_by(cancer_type) %>%
    summarise(
      total_genes_tested = n(),
      significant_associations = sum(p_value < 0.05),
      co_occurring = sum(p_value < 0.05 & odds_ratio > 1),
      mutually_exclusive = sum(p_value < 0.05 & odds_ratio < 1),
      .groups = 'drop'
    )
  
  print("Co-occurrence Analysis Summary:")
  print(summary_stats)
}

# --- Part 3 ---
# --- 3.1: 資料準備 - 整合生物標記物 ---
# TMB 和 MSI 來自基因檢測資料庫
biomarkers_tmb_msi <- report_biomarker_findings %>%
  pivot_wider(names_from = report_property, values_from = value) %>%
  rename(tmb_score = `TumorMutationBurdenScore`, msi_status = `Microsatellite status`) %>%
  mutate(
    tmb_score = as.numeric(str_remove(tmb_score, " Muts/Mb")),
    msi_group = case_when(
      msi_status == "MSI-High" ~ "MSI-High",
      msi_status == "MS-Stable" ~ "MSS",
      TRUE ~ "Unknown"
    )
  )

# PD-L1 來自臨床資料庫
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

# 將所有生物標記物與 MTAP 亞型資訊整合
analysis_fig3 <- master_df %>%
  select(hash_id, submitted_diagnosis, mtap_status, Cancer_type) %>%
  left_join(mtap_isoforms, by = "hash_id") %>%
  left_join(biomarkers_tmb_msi, by = "hash_id") %>%
  left_join(biomarkers_pd_l1, by = "hash_id") %>%
  # 建立MTAP分組
  mutate(
    mtap_group = case_when(
      mtap_status == "WT" ~ "MTAP Wild-type",
      !is.na(isoform_simplified) ~ "MTAP Altered",
      TRUE ~ "Unknown"
    ),
    # 根據Cancer_type編碼轉換癌別名稱
    cancer_name = case_when(
      Cancer_type == "01" ~ "Non-squamous NSCLC",
      Cancer_type == "02" ~ "Esophageal cancer",
      Cancer_type == "03" ~ "Gastric cancer",
      Cancer_type == "04" ~ "Gallbladder cancer",
      Cancer_type %in% c("05", "07", "08") ~ "Cholangiocarcinoma",
      Cancer_type == "06" ~ "Pancreatic cancer",
      TRUE ~ submitted_diagnosis
    )
  )

# --- 3.2: 繪製 Figure 3a - TMB 分佈比較 ---
create_tmb_plot <- function(cancer_filter = NULL) {
  plot_data <- analysis_fig3
  
  if (!is.null(cancer_filter)) {
    plot_data <- plot_data %>% filter(cancer_name == cancer_filter)
  }
  
  # 統計檢驗
  stat_test <- plot_data %>%
    filter(!is.na(tmb_score), !is.na(mtap_group)) %>%
    wilcox_test(tmb_score ~ mtap_group) %>%
    add_significance() %>%
    add_xy_position(x = "mtap_group")
  
  p <- ggplot(plot_data %>% filter(!is.na(tmb_score), !is.na(mtap_group)),
              aes(x = mtap_group, y = tmb_score, fill = mtap_group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif") +
    labs(
      title = ifelse(is.null(cancer_filter), "TMB Distribution by MTAP Status", 
                     paste("TMB Distribution by MTAP Status -", cancer_filter)),
      x = "MTAP Status",
      y = "TMB (Mutations/Mb)",
      fill = "MTAP Status"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("MTAP Wild-type" = "#66C2A5", "MTAP Altered" = "#FC8D62")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# 繪製各癌別的TMB分佈
cancer_types_fig3 <- c("Non-squamous NSCLC", "Pancreatic cancer", "Gastric cancer")
for (cancer_type in cancer_types_fig3) {
  print(create_tmb_plot(cancer_type))
}

# --- 3.3: 繪製 Figure 3b - PD-L1 表現分佈 ---
create_pdl1_plot <- function(cancer_filter = NULL) {
  plot_data <- analysis_fig3
  
  if (!is.null(cancer_filter)) {
    plot_data <- plot_data %>% filter(cancer_name == cancer_filter)
  }
  
  pd_l1_summary <- plot_data %>%
    filter(!is.na(pd_l1_group), !is.na(mtap_group)) %>%
    count(mtap_group, pd_l1_group) %>%
    group_by(mtap_group) %>%
    mutate(
      proportion = n / sum(n),
      total_n = sum(n)
    ) %>%
    ungroup()
  
  if (nrow(pd_l1_summary) > 0) {
    p <- ggplot(pd_l1_summary, aes(x = mtap_group, y = proportion, fill = pd_l1_group)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 3) +
      labs(
        title = ifelse(is.null(cancer_filter), "PD-L1 Expression by MTAP Status", 
                       paste("PD-L1 Expression by MTAP Status -", cancer_filter)),
        x = "MTAP Status",
        y = "Proportion",
        fill = "PD-L1 Expression"
      ) +
      scale_y_continuous(labels = scales::percent) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set2") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
}

# 繪製各癌別的PD-L1分佈
for (cancer_type in cancer_types_fig3) {
  p <- create_pdl1_plot(cancer_type)
  if (!is.null(p)) print(p)
}

# --- 3.4: 繪製 Figure 3c - MSI 狀態分佈 ---
create_msi_plot <- function(cancer_filter = NULL) {
  plot_data <- analysis_fig3
  
  if (!is.null(cancer_filter)) {
    plot_data <- plot_data %>% filter(cancer_name == cancer_filter)
  }
  
  msi_summary <- plot_data %>%
    filter(!is.na(msi_group), !is.na(mtap_group)) %>%
    count(mtap_group, msi_group) %>%
    group_by(mtap_group) %>%
    mutate(
      proportion = n / sum(n),
      total_n = sum(n)
    ) %>%
    ungroup()
  
  if (nrow(msi_summary) > 0) {
    p <- ggplot(msi_summary, aes(x = mtap_group, y = proportion, fill = msi_group)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 3) +
      labs(
        title = ifelse(is.null(cancer_filter), "MSI Status by MTAP Status", 
                       paste("MSI Status by MTAP Status -", cancer_filter)),
        x = "MTAP Status",
        y = "Proportion",
        fill = "MSI Status"
      ) +
      scale_y_continuous(labels = scales::percent) +
      theme_minimal() +
      scale_fill_manual(values = c("MSI-High" = "#E41A1C", "MSS" = "#377EB8")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
}

# 繪製各癌別的MSI分佈
for (cancer_type in cancer_types_fig3) {
  p <- create_msi_plot(cancer_type)
  if (!is.null(p)) print(p)
}

# --- 3.5: 統計摘要表 ---
biomarker_summary <- analysis_fig3 %>%
  filter(!is.na(mtap_group)) %>%
  group_by(cancer_name, mtap_group) %>%
  summarise(
    n = n(),
    tmb_median = median(tmb_score, na.rm = TRUE),
    tmb_iqr = IQR(tmb_score, na.rm = TRUE),
    pd_l1_high_pct = mean(pd_l1_group == "High (≥50%)", na.rm = TRUE) * 100,
    msi_high_pct = mean(msi_group == "MSI-High", na.rm = TRUE) * 100,
    .groups = 'drop'
  )

print("Biomarker Summary by Cancer Type and MTAP Status:")
print(biomarker_summary)

# --- 3.6: 相關性分析 ---
# 檢驗TMB與MTAP狀態的關聯性
correlation_results <- analysis_fig3 %>%
  filter(!is.na(tmb_score), !is.na(mtap_group)) %>%
  group_by(cancer_name) %>%
  group_modify(~ {
    if (length(unique(.x$mtap_group)) > 1) {
      test_result <- wilcox.test(tmb_score ~ mtap_group, data = .x)
      tibble(
        test = "Wilcoxon rank-sum test",
        p_value = test_result$p.value,
        n_wt = sum(.x$mtap_group == "MTAP Wild-type"),
        n_alt = sum(.x$mtap_group == "MTAP Altered"),
        median_tmb_wt = median(.x$tmb_score[.x$mtap_group == "MTAP Wild-type"], na.rm = TRUE),
        median_tmb_alt = median(.x$tmb_score[.x$mtap_group == "MTAP Altered"], na.rm = TRUE)
      )
    } else {
      tibble(
        test = "Insufficient data",
        p_value = NA,
        n_wt = NA,
        n_alt = NA,
        median_tmb_wt = NA,
        median_tmb_alt = NA
      )
    }
  })

print("TMB vs MTAP Status Correlation Analysis:")
print(correlation_results)