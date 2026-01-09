# Snippets Repo 規劃

## 目的

建立一個統一的 GitHub repository 來累積可重用的程式碼片段，取代散落在 Obsidian 筆記中難以維護的 code blocks。

## 目錄結構

```
snippets/
├── R/
│   ├── tidyverse/          # dplyr, ggplot2, tidyr 常用操作
│   ├── deseq2/             # RNA-seq 差異表達分析
│   ├── wgcna/              # 共表達網絡分析
│   ├── survival/           # 存活分析、CoxPH
│   ├── enrichment/         # GSEA, GO, ClusterProfiler
│   ├── visualization/      # pheatmap, 熱圖, 火山圖
│   └── stats/              # 統計檢定、ANOVA、post-hoc
│
├── python/
│   ├── ml/                 # sklearn, XGBoost, 模型訓練
│   ├── data-processing/    # pandas, encoding, SMOTE
│   ├── api/                # requests, API 呼叫
│   └── plotting/           # matplotlib, seaborn
│
├── bash/
│   ├── genomics/           # PLINK, SHAPEIT, IMPUTE5 指令
│   ├── ngs-pipeline/       # fastq處理, alignment, variant calling
│   └── server/             # SSH, SCP, job submission
│
├── sql/
│   └── queries/            # 常用查詢模板
│
├── docker/
│   └── templates/          # Dockerfile, docker-compose 範例
│
└── templates/
    ├── project-structure/  # 新專案資料夾結構
    └── rmarkdown/          # 報告模板
```

## 設計原則

1. **按語言 → 領域**：先依程式語言分類，再依應用場景細分
2. **扁平化**：最多兩層目錄，避免過度巢狀
3. **對應工作流**：涵蓋 NGS 分析、機器學習、統計、視覺化等日常任務

## 檔案命名規範

使用 kebab-case 描述性名稱：

```
R/deseq2/
├── basic-deseq2-workflow.R
├── lfc-shrinkage.R
└── batch-correction-with-svaseq.R

python/ml/
├── train-test-split.py
├── xgboost-hyperparameter-tuning.py
└── cross-validation-template.py
```

## Snippet 檔案格式建議

每個 snippet 檔案開頭加上簡短說明：

```r
# ---------------------------------------------
# Title: Basic DESeq2 Workflow
# Description: 從 count matrix 到差異表達結果
# Input: count matrix, sample metadata
# Output: DESeqDataSet, results table
# ---------------------------------------------
```

## 與 Obsidian 整合

在 Atlas 筆記中嵌入 GitHub 連結，保持筆記輕量：

```markdown
## DESeq2 差異表達分析

詳細程式碼：[GitHub](https://github.com/username/snippets/blob/main/R/deseq2/basic-deseq2-workflow.R)

重點步驟：
- 建立 DESeqDataSet
- 執行差異分析
- LFC shrinkage
```

## 下一步

- [ ] 建立 GitHub repository
- [ ] 初始化目錄結構
- [ ] 從現有 Obsidian 筆記遷移程式碼
- [ ] 設定 README.md 作為索引
