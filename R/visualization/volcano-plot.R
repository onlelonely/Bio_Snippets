# ---------------------------------------------
# Title: Volcano Plot for Co-occurrence Analysis
# Description: Generate volcano plots showing odds ratios vs p-values
# Input: Co-occurrence results with log2_or and log10_p
# Output: ggplot volcano plot object
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(ggplot2)
library(ggrepel)

#' Create volcano plot for co-occurrence analysis
plot_volcano <- function(cooccurrence_results, 
                         p_threshold = 0.05,
                         top_n_label = 10,
                         title = "Mutation Co-occurrence Analysis") {
  
  # Identify top genes to label
  top_cooccurring <- cooccurrence_results %>%
    filter(direction == "Co-occurring") %>%
    arrange(p_value) %>%
    head(top_n_label)
  
  top_exclusive <- cooccurrence_results %>%
    filter(direction == "Mutually exclusive") %>%
    arrange(p_value) %>%
    head(top_n_label)
  
  genes_to_label <- bind_rows(top_cooccurring, top_exclusive)$gene
  
  # Prepare plot data
  plot_data <- cooccurrence_results %>%
    mutate(
      label = ifelse(gene %in% genes_to_label, gene, ""),
      color_group = case_when(
        direction == "Co-occurring" ~ "Co-occurring (p < 0.05)",
        direction == "Mutually exclusive" ~ "Mutually exclusive (p < 0.05)",
        TRUE ~ "Not significant"
      )
    )
  
  ggplot(plot_data, aes(x = log2_or, y = log10_p, color = color_group)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray40") +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    scale_color_manual(
      values = c(
        "Co-occurring (p < 0.05)" = "firebrick",
        "Mutually exclusive (p < 0.05)" = "steelblue",
        "Not significant" = "gray60"
      )
    ) +
    labs(
      title = title,
      x = "log2(Odds Ratio)",
      y = "-log10(P-value)",
      color = "Association Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}
