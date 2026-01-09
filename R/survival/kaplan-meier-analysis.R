# ---------------------------------------------
# Title: Kaplan-Meier Survival Analysis
# Description: Fit KM model, log-rank test, and generate survival curves
# Input: Survival data with OS_time, OS_status
# Output: Survival object, KM fit, log-rank test results
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(survival)
library(survminer)

#' Fit Kaplan-Meier model and perform log-rank test
fit_km_analysis <- function(survival_data, group_var = "target_gene_status") {
  # Create survival object
  surv_object <- Surv(
    time = survival_data$OS_time,
    event = survival_data$OS_status
  )
  
  # Fit KM model by group
  formula <- as.formula(paste("surv_object ~", group_var))
  km_fit <- survfit(formula, data = survival_data)
  
  # Log-rank test
  log_rank <- survdiff(formula, data = survival_data)
  log_rank_p <- 1 - pchisq(log_rank$chisq, df = 1)
  
  list(
    surv_object = surv_object,
    km_fit = km_fit,
    log_rank = log_rank,
    log_rank_p = log_rank_p
  )
}

#' Generate Kaplan-Meier plot with survminer
plot_km_curves <- function(km_fit, survival_data, 
                           title = "Kaplan-Meier Survival Curves",
                           legend_title = "Group",
                           legend_labs = NULL,
                           palette = c("lightgray", "firebrick")) {
  ggsurvplot(
    km_fit,
    data = survival_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = palette,
    legend.title = legend_title,
    legend.labs = legend_labs,
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = title,
    ggtheme = theme_minimal()
  )
}
