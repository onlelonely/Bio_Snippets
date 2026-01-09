# ---------------------------------------------
# Title: Prepare Survival Analysis Data
# Description: Create survival time and status variables from clinical data
# Input: Master dataframe, death data, follow-up data
# Output: Dataframe with OS_time and OS_status columns
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(tidyverse)

#' Prepare survival analysis data
prepare_survival_data <- function(master_df, death_data, followup_data, columns, params) {
  key_col <- columns$key_column
  death_col <- columns$clinical$death_date
  followup_col <- columns$clinical$followup_date
  year_col <- columns$clinical$year_of_birth
  gender_col <- columns$clinical$gender
  
  days_per_month <- params$statistics$survival$days_per_month
  ref_year <- params$reference$year
  
  survival_df <- master_df %>%
    left_join(death_data %>% select(!!sym(key_col), !!sym(death_col)), by = key_col) %>%
    left_join(followup_data %>% select(!!sym(key_col), !!sym(followup_col)), by = key_col) %>%
    mutate(
      death_days = as.numeric(!!sym(death_col)),
      followup_days = as.numeric(!!sym(followup_col)),
      # OS_time in months
      OS_time = ifelse(!is.na(death_days), death_days, followup_days) / days_per_month,
      # OS_status: 1 = event (death), 0 = censored
      OS_status = ifelse(!is.na(death_days), 1, 0),
      # Demographics
      age_at_diagnosis = ref_year - as.numeric(!!sym(year_col)),
      gender_factor = factor(!!sym(gender_col), levels = c("1", "2"), labels = c("Male", "Female"))
    ) %>%
    filter(OS_time > 0, !is.na(OS_time))
  
  survival_df
}
