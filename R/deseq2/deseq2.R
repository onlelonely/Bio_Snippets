# ---------------------------------------------
# Title: DESeq2
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/DESeq2.md
# ---------------------------------------------

### Execute differential expression analysis using DESeq2
mrna_dds <-
  DESeqDataSetFromMatrix(round(mrna_df),
                         colData = mrna_meta,
                         design = ~ Condition)

mrna_dds$Condition <- relevel(mrna_dds$Condition, ref = "Normal")
mrna_dds <- DESeq(mrna_dds)
vsd <- varianceStabilizingTransformation(mrna_dds, blind = FALSE)
head(assay(vsd))
hist(assay(vsd))
resultsNames(mrna_dds)
# Plot Dispersions:
plotDispEsts(mrna_dds, main = "Dispersion plot")
# MA plot
mrna_res <- results(mrna_dds, name = "Condition_Tumor_vs_Normal")
plotMA(mrna_res)

#Statistics summary
mrna_res_df <- as.data.frame(mrna_res)
mrnaTable <- mrna_res_df
mrnaTable$Gene_id <- rownames(mrnaTable)
summary(mrna_res)						 

#Volcuno plot
mrna_upreg <- get_upregulated(mrna_res)
mrna_downreg <- get_downregulated(mrna_res)
mrna_counts <- counts(mrna_dds, normalized = T)
mrna_upreg$Gene_id <- rownames(mrna_upreg)
mrna_downreg$Gene_id <- rownames(mrna_downreg)
mrna_res_df$Gene_id <- rownames(mrna_res_df)		

par(
  mar = c(5, 5, 5, 5),
  cex = 1.0,
  cex.main = 1.4,
  cex.axis = 1.4,
  cex.lab = 1.4
)
with(
  mrna_res_df,
  plot(
    log2FoldChange,
    -log10(padj),
    pch = 20,
    main = "Volcano plot",
    cex = 1.0,
    xlab = bquote( ~ Log[2] ~ fold ~ change),
    ylab = bquote( ~ -log[10] ~ P ~ value)
  )
)

with(
  subset(mrna_res_df, padj < 0.05 &
           abs(log2FoldChange) > 2),
  points(
    log2FoldChange,
    -log10(padj),
    pch = 20,
    col = "red",
    cex = 0.5
  )
)

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v = 0,
       col = "black",
       lty = 3,
       lwd = 1.0)
abline(v = -2,
       col = "black",
       lty = 4,
       lwd = 2.0)
abline(v = 2,
       col = "black",
       lty = 4,
       lwd = 2.0)
abline(
  h = -log10(max(mrna_res_df$pvalue[mrna_res_df$padj < 0.05], na.rm = TRUE)),
  col = "black",
  lty = 4,
  lwd = 2.0
)

rownames(mrna_upreg) <- gsub("\\..*", "", rownames(mrna_upreg))
rownames(mrna_downreg) <- gsub("\\..*", "", rownames(mrna_downreg))

Genes <- read.csv("geneList.csv", header=FALSE)
library(Biomart)
Ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
results = getBM(
    attributes = c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'), 
    filters = 'hgnc_symbol',
    values = Genes,
    mart = ensembl
)
library(dplyr)
mrna_upreg_filtered <- mrna_upreg[rownames(mrna_upreg) %in% results$ensembl_gene_id,]
mrna_downreg_filtered <- mrna_upreg[rownames(mrna_downreg) %in% results$ensembl_gene_id,]
sig_genes <- rbind(mrna_upreg_filtered, mrna_downreg_filtered)