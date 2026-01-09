# ---------------------------------------------
# Title: Cox Proportional Hazards Regression
# Description: Univariate and multivariate Cox models with forest plots
# Input: Survival data with covariates
# Output: Cox model results, hazard ratios, forest plot
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(survival)
library(broom)
library(ggplot2)

#' Fit multivariate Cox model
fit_cox_model <- function(surv_object, survival_data, covariates) {
  # Build formula
  formula <- as.formula(paste("surv_object ~", paste(covariates, collapse = " + ")))
  
  # Fit model
  cox_fit <- coxph(formula, data = survival_data)
  
  # Extract hazard ratios
  hr_table <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(
      Term = term,
      HR = estimate,
      lower_ci = conf.low,
      upper_ci = conf.high,
      p_value = p.value
    )
  
  list(
    model = cox_fit,
    hr_table = hr_table
  )
}

#' Test proportional hazards assumption
test_ph_assumption <- function(cox_fit) {
  cox.zph(cox_fit)
}

#' Create forest plot of hazard ratios
plot_forest <- function(hr_table, title = "Cox Regression: Hazard Ratios") {
  hr_table %>%
    mutate(significant = p_value < 0.05) %>%
    ggplot(aes(x = HR, y = reorder(Term, HR), color = significant)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    scale_color_manual(
      values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
      labels = c("Not significant", "p < 0.05")
    ) +
    scale_x_log10() +
    labs(title = title, x = "Hazard Ratio (log scale)", y = "", color = "Significance") +
    theme_minimal() +
    theme(legend.position = "bottom")
}
