# ---------------------------------------------
# Title: Exploratory Data Analysis Plots
# Description: Common EDA visualizations - histograms, boxplots, bar charts
# Input: Dataframe with variables to visualize
# Output: ggplot objects
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(ggplot2)
library(ggpubr)

#' Age distribution histogram with median line
plot_age_distribution <- function(df, age_col = "age") {
  median_age <- median(dfage_col, na.rm = TRUE)
  
  ggplot(df, aes(x = .dataage_col)) +
    geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
    geom_vline(aes(xintercept = median_age), color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Age Distribution",
      subtitle = paste("Median:", median_age, "years"),
      x = "Age (years)",
      y = "Count"
    ) +
    theme_minimal()
}

#' BMI distribution with category thresholds
plot_bmi_distribution <- function(df, bmi_col = "bmi") {
  ggplot(df, aes(x = .databmi_col)) +
    geom_histogram(binwidth = 2, fill = "coral", color = "white") +
    geom_vline(xintercept = 25, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
    labs(
      title = "BMI Distribution",
      subtitle = "Blue: Normal/Overweight (25), Red: Overweight/Obese (30)",
      x = "BMI",
      y = "Count"
    ) +
    theme_minimal()
}

#' Horizontal bar chart with percentages
plot_category_bars <- function(df, category_col, title = "Distribution") {
  summary_df <- df %>%
    count(.datacategory_col) %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))
  
  ggplot(summary_df, aes(x = reorder(.datacategory_col, n), y = n, fill = .datacategory_col)) +
    geom_col() +
    geom_text(aes(label = paste0(n, " (", pct, "%)")), hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(title = title, x = "", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
}

#' Grouped boxplot comparison
plot_grouped_boxplot <- function(df, x_col, y_col, fill_col = NULL,
                                  title = "Comparison", x_lab = "", y_lab = "") {
  p <- ggplot(df, aes(x = .datax_col, y = .datay_col))
  
  if (!is.null(fill_col)) {
    p <- p + aes(fill = .datafill_col)
  }
  
  p +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, width = 0.2) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_minimal()
}
