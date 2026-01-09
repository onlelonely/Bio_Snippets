# ---------------------------------------------
# Title: YAML Configuration Loader
# Description: Load and parse YAML config files for R analysis pipelines
# Input: YAML files in config directory
# Output: List containing configuration objects
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
})

#' Load main project configuration
load_config <- function(config_dir = "config") {
  config_file <- file.path(config_dir, "config.yaml")
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }
  yaml.load_file(config_file)
}

#' Load all configuration files at once
load_all_config <- function(config_dir = "config") {
  list(
    config = yaml.load_file(file.path(config_dir, "config.yaml")),
    columns = yaml.load_file(file.path(config_dir, "column_mappings.yaml")),
    cancer_types = yaml.load_file(file.path(config_dir, "cancer_types.yaml")),
    analysis_params = yaml.load_file(file.path(config_dir, "analysis_params.yaml"))
  )
}

#' Get nested parameter value
get_param <- function(params, ...) {
  param_path <- list(...)
  result <- params
  for (key in param_path) {
    if (is.null(resultkey)) {
      stop(paste("Parameter not found:", paste(param_path, collapse = " -> ")))
    }
    result <- resultkey
  }
  result
}
