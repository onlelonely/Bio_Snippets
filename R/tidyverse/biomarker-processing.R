# ---------------------------------------------
# Title: Biomarker Data Processing
# Description: Process TMB, MSI, and PD-L1 biomarker data
# Input: Raw biomarker dataframes
# Output: Processed biomarker data with categories
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(tidyverse)

#' Process TMB and MSI biomarker findings (long to wide)
process_biomarker_findings <- function(biomarkers, columns) {
  tmb_name <- columns$biomarker_findings$tmb_property_name
  msi_name <- columns$biomarker_findings$msi_property_name
  prop_col <- columns$biomarker_findings$report_property
  value_col <- columns$biomarker_findings$value
  key_col <- columns$key_column
  
  biomarkers_wide <- biomarkers %>%
    pivot_wider(
      id_cols = !!sym(key_col),
      names_from = !!sym(prop_col),
      values_from = !!sym(value_col)
    )
  
  # Rename and process TMB
  if (tmb_name %in% colnames(biomarkers_wide)) {
    biomarkers_wide <- biomarkers_wide %>%
      rename(tmb_score = !!sym(tmb_name)) %>%
      mutate(tmb_score = as.numeric(str_remove(tmb_score, " Muts/Mb")))
  }
  
  # Rename MSI
  if (msi_name %in% colnames(biomarkers_wide)) {
    biomarkers_wide <- biomarkers_wide %>%
      rename(msi_status = !!sym(msi_name))
  }
  
  biomarkers_wide
}

#' Categorize PD-L1 expression levels
process_pdl1_data <- function(biomarker_data, columns, params) {
  key_col <- columns$key_column
  pdl1_percent_col <- columns$biomarker$pdl1_percent
  high_threshold <- params$biomarkers$pdl1$high_threshold
  low_threshold <- params$biomarkers$pdl1$low_threshold
  
  biomarker_data %>%
    mutate(
      pd_l1_value = as.numeric(!!sym(pdl1_percent_col)),
      pd_l1_group = case_when(
        pd_l1_value >= high_threshold ~ sprintf("High (≥%d%%)", high_threshold),
        pd_l1_value >= low_threshold ~ sprintf("Low (%d-%d%%)", low_threshold, high_threshold - 1),
        pd_l1_value < low_threshold ~ sprintf("Negative (<%d%%)", low_threshold),
        TRUE ~ "Unknown"
      )
    ) %>%
    select(!!sym(key_col), pd_l1_group, pd_l1_value)
}
