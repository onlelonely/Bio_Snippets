# ---------------------------------------------
# Title: YAML Configuration Loader
# Description: Load multiple YAML config files for Python pipelines
# Input: Config directory with YAML files
# Output: Dictionary with configuration objects
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

import yaml
import os


def load_all_config(config_dir="config"):
    """Load all YAML configuration files."""
    config = {}
    
    config_files = {
        'config': 'config.yaml',
        'columns': 'column_mappings.yaml',
        'cancer_types': 'cancer_types.yaml',
        'params': 'analysis_params.yaml'
    }
    
    for key, filename in config_files.items():
        filepath = os.path.join(config_dir, filename)
        with open(filepath, 'r') as f:
            config[key] = yaml.safe_load(f)
    
    return config


def get_nested_param(config, *keys):
    """Get nested parameter value from config dict."""
    result = config
    for key in keys:
        if key not in result:
            raise KeyError(f"Key not found: {' -> '.join(keys)}")
        result = result[key]
    return result
