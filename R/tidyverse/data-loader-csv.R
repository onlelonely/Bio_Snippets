# ---------------------------------------------
# Title: CSV Data Loader with Config
# Description: Load multiple CSV files using configuration-based paths
# Input: Config object with file paths, column mappings
# Output: List of data frames
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

#' Load clinical data files from config
load_clinical_data <- function(config, columns) {
  data_dir <- file.path(config$paths$data_dir, config$paths$clinical_subdir)
  
  if (!dir.exists(data_dir)) {
    stop(paste("Data directory not found:", data_dir))
  }
  
  list(
    demographics = read_csv(file.path(data_dir, config$files$clinical$demographics), show_col_types = FALSE),
    cancer_char = read_csv(file.path(data_dir, config$files$clinical$cancer_char), show_col_types = FALSE),
    biomarker = read_csv(file.path(data_dir, config$files$clinical$biomarker), show_col_types = FALSE),
    death = read_csv(file.path(data_dir, config$files$clinical$death), show_col_types = FALSE),
    followup = read_csv(file.path(data_dir, config$files$clinical$followup), show_col_types = FALSE)
  )
}

#' Load genomic data files from config
load_genomic_data <- function(config, columns) {
  data_dir <- file.path(config$paths$data_dir, config$paths$genomic_subdir)
  
  list(
    report = read_csv(file.path(data_dir, config$files$genomic$report), show_col_types = FALSE),
    biomarkers = read_csv(file.path(data_dir, config$files$genomic$biomarkers), show_col_types = FALSE),
    snv = read_csv(file.path(data_dir, config$files$genomic$snv), show_col_types = FALSE),
    cna = read_csv(file.path(data_dir, config$files$genomic$cna), show_col_types = FALSE),
    rearrangement = read_csv(file.path(data_dir, config$files$genomic$rearrangement), show_col_types = FALSE)
  )
}
