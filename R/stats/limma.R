# ---------------------------------------------
# Title: Limma
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/Limma.md
# ---------------------------------------------

library(edgeR)
library(limma)
library(Biomart)

geneList <- read.csv("./Data/DryLab_P6/P6_geneList.txt")
Ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

results = getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'), 
  filters = 'hgnc_symbol',
  values = geneList,
  mart = ensembl
)
# Apply gene list
cleaned_rownames <- gsub("\\..*", "", rownames(mrna_df))
mrna_df_filtered <- mrna_df[cleaned_rownames %in% results$ensembl_gene_id, ]
# Filter lowly expressed genes
filter_condition <- rowSums(mrna_df_filtered >= 1) > (ncol(mrna_df_filtered) / 2)
mrna_df_filtered <- mrna_df_filtered[filter_condition, ]

# Create a design matrix
design <- model.matrix(~ 0 + Condition, data = mrna_meta)
colnames(design) <- levels(mrna_meta$Condition)

# Create a DGEList object
y <- DGEList(counts=mrna_df_filtered)
# Calculate normalization factors
y <- calcNormFactors(y)
# Apply voom transformation
v <- voom(y, design, plot=TRUE)
# Fit the linear model
fit <- lmFit(v, design)

# Proceed with contrasts and other limma steps as before
cont.matrix <- makeContrasts(Tumor-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get the results table
DEG <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# Filter genes by some criteria
sig.genes <- DEG[DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > 1.5,]