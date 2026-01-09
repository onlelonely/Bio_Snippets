# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

A unified collection of reusable code snippets for bioinformatics and data science workflows, replacing scattered code blocks in Obsidian notes.

## Architecture

**Directory structure**: Language → Domain (max 2 levels)
- `R/` - tidyverse, DESeq2, WGCNA, survival analysis, enrichment, visualization, statistics
- `python/` - ML, data processing, API, plotting
- `bash/` - genomics (PLINK, SHAPEIT, IMPUTE5), NGS pipelines, server commands
- `sql/` - query templates
- `docker/` - Dockerfile and docker-compose templates
- `templates/` - project structures, RMarkdown templates

## Conventions

**File naming**: kebab-case (e.g., `basic-deseq2-workflow.R`, `xgboost-hyperparameter-tuning.py`)

**Snippet header format**:
```
# ---------------------------------------------
# Title: [Descriptive title]
# Description: [What it does]
# Input: [Expected inputs]
# Output: [What it produces]
# ---------------------------------------------
```
