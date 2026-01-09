# Bio_Snippets

A collection of reusable code snippets for bioinformatics and data science workflows.

## Repository Structure

```
Bio_Snippets/
├── R/
│   ├── deseq2/          # Differential expression analysis
│   ├── enrichment/      # GO, GSEA, pathway analysis
│   ├── stats/           # Statistical methods, clustering metrics
│   ├── survival/        # Kaplan-Meier, Cox regression
│   ├── tcga/            # TCGA data analysis workflows
│   ├── tidyverse/       # Data wrangling, config loading
│   ├── visualization/   # Plots, heatmaps, volcano plots
│   └── wgcna/           # Weighted gene co-expression network
├── python/
│   ├── api/             # API clients (AlphaFold, etc.)
│   ├── data-processing/ # Encoding, data wrangling
│   ├── drug-screening/  # Virtual screening pipeline
│   └── ml/              # Machine learning models
├── bash/
│   ├── genomics/        # Genomics tools
│   ├── ngs-pipeline/    # NGS analysis scripts
│   ├── server/          # Server utilities
│   └── setup/           # Environment setup scripts
├── docker/
│   └── templates/       # Docker compose templates
├── templates/
│   └── rmarkdown/       # Analysis notebook templates
└── sql/                 # SQL query templates
```

## Highlights

### R
- **Survival Analysis**: Kaplan-Meier curves, Cox proportional hazards, forest plots
- **Enrichment**: clusterProfiler, GSEA, ssGSEA, GO analysis
- **TCGA Workflows**: Data loading, DESeq2 integration, survival analysis
- **Visualization**: Volcano plots, heatmaps, consensus clustering

### Python
- **Drug Screening Pipeline**: 14-module virtual screening workflow
  - Receptor preparation, ligand validation, parallel docking
  - MD simulation, results analysis, report generation
- **ML Models**: SMOTE, target encoding, species-aware attention networks
- **Data Processing**: YAML config loaders, simulated data generators

### Bash
- **NGS Pipelines**: HLA typing, T-cell analysis, VDJ recombination
- **Setup Scripts**: ColabFold, AI-Scientist, environment configuration

## Snippet Format

Each snippet follows a standard header format:

```r
# ---------------------------------------------
# Title: Descriptive Title
# Description: What it does
# Input: Expected inputs
# Output: What it produces
# ---------------------------------------------
```

## Usage

Browse the folders and copy the snippets you need. Most snippets are self-contained functions that can be adapted to your workflow.

## License

MIT License
