# ---------------------------------------------
# Title: Variant Data Standardization
# Description: Combine SNV, CNA, rearrangement data into unified format
# Input: Separate variant dataframes
# Output: Combined standardized dataframe
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

library(tidyverse)

#' Standardize variant data from multiple sources
standardize_variants <- function(snv, cna, rearrangements, columns) {
  key_col <- columns$key_column
  gene_col <- columns$genomic$gene
  protein_col <- columns$genomic$protein_effect
  type_col <- columns$genomic$type
  targeted_col <- columns$genomic$targeted_gene
  other_col <- columns$genomic$other_gene
  
  # Standardize short variants (SNV)
  sv_std <- snv %>%
    transmute(
      !!sym(key_col) := !!sym(key_col),
      gene = !!sym(gene_col),
      alteration_type = "SV",
      alteration_detail = !!sym(protein_col)
    )
  
  # Standardize copy number alterations (CNA)
  cna_std <- cna %>%
    transmute(
      !!sym(key_col) := !!sym(key_col),
      gene = !!sym(gene_col),
      alteration_type = "CNA",
      alteration_detail = !!sym(type_col)
    )
  
  # Standardize rearrangements (RE)
  re_std <- rearrangements %>%
    rename(gene = !!sym(targeted_col)) %>%
    transmute(
      !!sym(key_col) := !!sym(key_col),
      gene = gene,
      alteration_type = "RE",
      alteration_detail = paste("fusion with", !!sym(other_col))
    )
  
  bind_rows(sv_std, cna_std, re_std)
}

#' Create gene mutation matrix (wide format)
create_gene_matrix <- function(variants_combined, patient_ids, columns) {
  key_col <- columns$key_column
  gene_col <- "gene"
  
  variants_combined %>%
    filter(!!sym(key_col) %in% patient_ids) %>%
    select(!!sym(key_col), !!sym(gene_col)) %>%
    filter(!is.na(!!sym(gene_col))) %>%
    distinct() %>%
    mutate(mutated = 1) %>%
    pivot_wider(names_from = !!sym(gene_col), values_from = mutated, values_fill = 0)
}
