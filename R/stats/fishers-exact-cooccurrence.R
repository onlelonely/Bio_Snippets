# ---------------------------------------------
# Title: Fisher's Exact Test for Co-occurrence Analysis
# Description: Analyze gene co-occurrence patterns using Fisher's exact test
# Input: Gene mutation matrix (patients x genes)
# Output: Dataframe with odds ratios, p-values, and significance
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(tidyverse)

#' Perform co-occurrence analysis using Fisher's exact test
perform_cooccurrence_analysis <- function(gene_matrix, target_gene, key_col, 
                                          min_mutations = 3, p_threshold = 0.05) {
  results <- list()
  other_genes <- setdiff(colnames(gene_matrix), c(key_col, target_gene))
  
  for (gene_name in other_genes) {
    # Skip genes with too few mutations
    gene_mut_count <- sum(gene_matrixgene_name)
    if (gene_mut_count < min_mutations) next
    
    # Create 2x2 contingency table
    contingency_table <- table(
      Target = gene_matrixtarget_gene,
      Other = gene_matrixgene_name
    )
    
    if (nrow(contingency_table) == 2 && ncol(contingency_table) == 2) {
      test_result <- fisher.test(contingency_table)
      
      # Calculate counts
      both <- sum(gene_matrixtarget_gene == 1 & gene_matrixgene_name == 1)
      target_only <- sum(gene_matrixtarget_gene == 1 & gene_matrixgene_name == 0)
      gene_only <- sum(gene_matrixtarget_gene == 0 & gene_matrixgene_name == 1)
      neither <- sum(gene_matrixtarget_gene == 0 & gene_matrixgene_name == 0)
      
      resultsgene_name <- tibble(
        gene = gene_name,
        p_value = test_result$p.value,
        odds_ratio = as.numeric(test_result$estimate),
        ci_low = test_result$conf.int[1],
        ci_high = test_result$conf.int[2],
        both_mutated = both,
        target_only = target_only,
        gene_only = gene_only,
        neither = neither,
        total_mutations = gene_mut_count
      )
    }
  }
  
  bind_rows(results) %>%
    mutate(
      log10_p = -log10(p_value),
      log2_or = log2(odds_ratio),
      significant = p_value < p_threshold,
      direction = case_when(
        !significant ~ "Not significant",
        odds_ratio > 1 ~ "Co-occurring",
        odds_ratio < 1 ~ "Mutually exclusive",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(p_value)
}

#' Display contingency table with Fisher's test results
show_contingency <- function(gene_matrix, gene_name, target_gene) {
  ct <- table(
    Target = gene_matrixtarget_gene,
    Other = gene_matrixgene_name
  )
  
  fisher_result <- fisher.test(ct)
  
  list(
    table = addmargins(ct),
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    conf_int = fisher_result$conf.int
  )
}
