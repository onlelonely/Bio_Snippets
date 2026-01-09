# ---------------------------------------------
# Title: TCGA analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/TCGA analysis.md
# ---------------------------------------------

library(TCGAbiolinks)
library(dplyr)
library(stringr) 
library(ggplot2)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(ggrepel)

### Previous code continues from mutation grouping...
# Assuming we have mrna_df_plain and mrna_meta from previous steps

### Load functions for identifying upregulated and downregulated genes
get_upregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange >= 1)],
                   rownames(df)[which(df$pvalue <= 0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key), ])
  return(results)
}

get_downregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange <= -1)],
                   rownames(df)[which(df$pvalue <= 0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key), ])
  return(results)
}

### Prepare data for DESeq2
# Remove genes with low expression (filter out genes with less than 10 reads across all samples)
keep <- rowSums(mrna_df_plain) >= 10
mrna_df_filtered <- mrna_df_plain[keep, ]

# Ensure count data is integer
mrna_df_filtered <- round(mrna_df_filtered)

# Match sample names between expression data and metadata
common_samples <- intersect(colnames(mrna_df_filtered), mrna_meta$cases)
mrna_df_filtered <- mrna_df_filtered[, common_samples]
mrna_meta <- mrna_meta[mrna_meta$cases %in% common_samples, ]

# Reorder metadata to match expression data
mrna_meta <- mrna_meta[match(colnames(mrna_df_filtered), mrna_meta$cases), ]

print(paste("Final sample count:", ncol(mrna_df_filtered)))
print("TP53 mutation status distribution:")
table(mrna_meta$Condition)

### Create DESeq2 dataset
# Set TP53_Wild_Type as reference level
mrna_meta$Condition <- factor(mrna_meta$Condition, levels = c("TP53_Wild_Type", "TP53_Mutated"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mrna_df_filtered,
                              colData = mrna_meta,
                              design = ~ Condition)

# Set reference level
dds$Condition <- relevel(dds$Condition, ref = "TP53_Wild_Type")

### Run DESeq2 analysis
print("Running DESeq2 analysis...")
dds <- DESeq(dds)

# Get results (TP53_Mutated vs TP53_Wild_Type)
res <- results(dds, contrast = c("Condition", "TP53_Mutated", "TP53_Wild_Type"))
res <- res[order(res$pvalue), ]

print("DESeq2 analysis completed!")
print(summary(res))

### Check TP53 expression specifically
# Find TP53 gene (might be listed under different gene name formats)
tp53_patterns <- c("TP53", "ENSG00000141510") # ENSEMBL ID for TP53
tp53_row <- NULL

for(pattern in tp53_patterns) {
  matches <- grep(pattern, rownames(res), ignore.case = TRUE)
  if(length(matches) > 0) {
    tp53_row <- matches[1]
    break
  }
}

if(!is.null(tp53_row)) {
  tp53_result <- res[tp53_row, ]
  print("=== TP53 Expression Analysis ===")
  print(paste("Gene ID:", rownames(res)[tp53_row]))
  print(paste("Log2 Fold Change:", round(tp53_result$log2FoldChange, 3)))
  print(paste("P-value:", format(tp53_result$pvalue, scientific = TRUE)))
  print(paste("Adjusted P-value:", format(tp53_result$padj, scientific = TRUE)))
  
  if(tp53_result$log2FoldChange < 0) {
    print("TP53 is DOWNREGULATED in TP53-mutated samples")
  } else {
    print("TP53 is UPREGULATED in TP53-mutated samples")
  }
} else {
  print("Warning: TP53 gene not found in the dataset")
  print("Available gene IDs (first 20):")
  print(head(rownames(res), 20))
}

### Extract significant differentially expressed genes
res_sig <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
print(paste("Number of significantly differentially expressed genes:", nrow(res_sig)))

# Get upregulated and downregulated genes
res_df <- as.data.frame(res_sig)
up_genes <- get_upregulated(res_df)
down_genes <- get_downregulated(res_df)

print(paste("Upregulated genes (log2FC >= 1, p < 0.05):", nrow(up_genes)))
print(paste("Downregulated genes (log2FC <= -1, p < 0.05):", nrow(down_genes)))

### Visualization 1: Volcano plot
res_df_all <- as.data.frame(res)
res_df_all$gene <- rownames(res_df_all)
res_df_all$significant <- ifelse(res_df_all$padj < 0.05 & !is.na(res_df_all$padj), "Yes", "No")
res_df_all$regulation <- ifelse(res_df_all$log2FoldChange > 1 & res_df_all$padj < 0.05, "Up",
                               ifelse(res_df_all$log2FoldChange < -1 & res_df_all$padj < 0.05, "Down", "NS"))

# Highlight TP53 if found
if(!is.null(tp53_row)) {
  res_df_all$highlight <- ifelse(rownames(res_df_all) == rownames(res)[tp53_row], "TP53", "Other")
} else {
  res_df_all$highlight <- "Other"
}

volcano_plot <- ggplot(res_df_all, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = regulation), alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Volcano Plot: TP53 Mutated vs TP53 Wild Type",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal()

# Add TP53 label if found
if(!is.null(tp53_row)) {
  tp53_data <- res_df_all[rownames(res_df_all) == rownames(res)[tp53_row], ]
  volcano_plot <- volcano_plot +
    geom_point(data = tp53_data, aes(x = log2FoldChange, y = -log10(pvalue)), 
               color = "orange", size = 3) +
    geom_text_repel(data = tp53_data, aes(label = "TP53"), 
                    color = "orange", fontface = "bold", size = 4)
}

print(volcano_plot)

### Visualization 2: Top differentially expressed genes heatmap
if(nrow(res_sig) > 0) {
  # Get top 50 most significant genes
  top_genes <- head(res_sig, 50)
  
  # Get normalized counts for heatmap
  vsd <- vst(dds, blind = FALSE)
  top_genes_matrix <- assay(vsd)[rownames(top_genes), ]
  
  # Create annotation for samples
  annotation_col <- data.frame(
    TP53_Status = mrna_meta$Condition
  )
  rownames(annotation_col) <- mrna_meta$cases
  
  # Create heatmap
  pheatmap(top_genes_matrix,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = FALSE,
           main = "Top 50 Differentially Expressed Genes\n(TP53 Mutated vs Wild Type)")
}

### Save results
write.csv(as.data.frame(res), "TP53_mutation_DESeq2_results.csv")
write.csv(up_genes, "TP53_mutation_upregulated_genes.csv")
write.csv(down_genes, "TP53_mutation_downregulated_genes.csv")

print("Analysis completed! Results saved to CSV files.")