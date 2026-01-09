# ---------------------------------------------
# Title: Amgen 報告Draft
# Description: From: Source/3. Efforts/Res_Amgen_MTAP_KRAS/_Archives/Amgen 報告Draft.md
# ---------------------------------------------

library(Tidyverse) 
library(survival)    
library(survminer)  
library(ggpubr)
library(purrr)
library(rstatix)  
library(ggrepel)  

# --- 載入資料 ---
# 統一由此修改載入的CSV檔名
# 臨床資料庫 (癌症主題資料檔)
patient_demographic <- read_csv("simulated_Clinical_Patient_Demographic.csv")
cancer_characteristic_new <- read_csv("simulated_Clinical_Cancer_Characteristic_Newly_Diagnosis.csv")
access_critical_biomarker <- read_csv("simulated_Clinical_Access_Cancer_Characteristic_Critical_Biomarker.csv")
death_data <- read_csv("simulated_Clinical_Death_Data.csv")

# 基因檢測資料庫 (癌症基因檢測檔)
report_patient <- read_csv("simulated_Genomic_Report_Patient.csv")
report_biomarker_findings <- read_csv("simulated_Genomic_Biomarker_Findings.csv")
variant_short <- read_csv("simulated_Genomic_Variant_Short-Variants.csv")
variant_cna <- read_csv("simulated_Genomic_Variant_Copy-Number-Alterations.csv")
variant_rearrangements <- read_csv("simulated_Genomic_Variant_Rearrangements.csv")

# --- 組合資料 ---
sv_standardized <- variant_short %>%
  transmute(
    hash_id,
    gene,
    alteration_type = "SV",
    alteration_detail = protein_effect 
  )

cna_standardized <- variant_cna %>%
  transmute(
    hash_id,
    gene,
    alteration_type = "CNA",
    alteration_detail = type 
  )

re_standardized <- variant_rearrangements %>%
  rename(gene = targeted_gene) %>%
  transmute(
    hash_id,
    gene,
    alteration_type = "RE",
    alteration_detail = paste("fusion with", othe_gene)
  )

variants_combined <- bind_rows(
  sv_standardized,
  cna_standardized,
  re_standardized
)

# 將基因檢測的病患資訊與臨床人口學資訊合併
master_df <- report_patient %>%
  left_join(patient_demographic, by = "hash_id") %>%
  left_join(cancer_characteristic_new, by = "hash_id")

# 建立 MTAP 狀態欄位
mtap_mutated_patients <- variants_combined %>%
  filter(gene == "MTAP") %>%
  distinct(hash_id)

# Grouping MTAP有無變異
master_df <- master_df %>%
  mutate(mtap_status = ifelse(hash_id %in% mtap_mutated_patients$hash_id, "Mutated", "WT"))