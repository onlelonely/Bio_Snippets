# ---------------------------------------------
# Title: GSEA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/GSEA.md
# ---------------------------------------------

library(pd.clariom.s.rat)
library(oligo)
library(affycoretools)
library(ggplot2)
library(ClusterProfiler)
library(GSVA)
library(msigdbr)
library(org.Rn.eg.db)
library(enrichplot)
library(pheatmap)
library(Biomart)

# 這段目前不能在docker跑，只能在local Rstudio跑
#celfiles <- list.files("../../permanent_storage/NTHU_Yeh_1/data/", 
#full = TRUE, pattern="*.CEL")
#celData <- read.celfiles(celfiles)
#eset <- oligo::RMA(celData, normalize=TRUE, background=TRUE)
#all.eset <- annotateEset(eset, annotation(eset))
#tc_rows <- grep("^TC", rownames(all.eset), value = TRUE)
#tc_eset <- all.eset[tc_rows, ]
#raw_expr <- exprs(tc_eset)
#probe_reference <- as.data.frame(fData(all.eset))

load("../../permanent_storage/NTHU_Yeh_1/data/post_RMA.Rdata")
eset_final <- raw_expr
indices <- match(rownames(eset_final), probe_reference$PROBEID)
symbols <- probe_reference$SYMBOL[indices]
# Add numbers to duplicate gene symbols
dup_symbols <- duplicated(symbols) | duplicated(symbols, fromLast=TRUE)
if(any(dup_symbols)) {
  for(sym in unique(symbols[dup_symbols])) {
    dups <- which(symbols == sym)
    symbols[dups] <- paste0(sym, "_", seq_along(dups))
  }
}
rownames(eset_final) <- symbols
colnames(eset_final) <- c("Control_1", "Control_2", "Control_3", "FUS_1", "FUS_2", "FUS_3")

# Remove rows that start with "---"
eset_final <- eset_final[!grepl("^---", rownames(eset_final)),]
# Find base gene names that have duplicates (ending in _number)
base_genes <- gsub("_\\d+$", "", rownames(eset_final))
dup_bases <- unique(base_genes[duplicated(base_genes)])
# For each duplicated base gene name
for(gene in dup_bases) {
  # Find all rows with this base gene
  dup_rows <- grep(paste0("^", gene, "_\\d+$"), rownames(eset_final))
  
  if(length(dup_rows) > 1) {
    # Calculate mean expression across duplicate entries
    mean_expr <- colMeans(eset_final[dup_rows,, drop=FALSE])
    
    # Remove duplicate rows
    eset_final <- eset_final[-dup_rows,]
    
    # Add back averaged row with base gene name
    eset_final <- rbind(eset_final, mean_expr)
    rownames(eset_final)[nrow(eset_final)] <- gene
  }
}

# Get significant genes
sig_genes <- results$ID[results$Significant]

# Convert gene symbols to ENTREZ IDs using org.Rn.eg.db (since this is rat data)
entrez_ids <- bitr(sig_genes, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID", 
                  OrgDb = "org.Rn.eg.db")

# Try mapping with biomaRt as an alternative approach
Ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl")

# Map using biomaRt
mart_ids <- getBM(attributes=c("external_gene_name", "entrezgene_id"),
                 filters="external_gene_name",
                 values=sig_genes,
                 mart=ensembl)

# Clean up biomaRt results
mart_ids <- mart_ids[!is.na(mart_ids$entrezgene_id),]
mart_ids <- unique(mart_ids)

# Compare mapping rates
cat("Mapping rate comparison:\n")
cat("org.Rn.eg.db mapped:", nrow(entrez_ids), "out of", length(sig_genes), "genes (",
    round(nrow(entrez_ids)/length(sig_genes)*100, 1), "%)\n")
cat("biomaRt mapped:", nrow(mart_ids), "out of", length(sig_genes), "genes (",
    round(nrow(mart_ids)/length(sig_genes)*100, 1), "%)\n")

# Compare overlapping and unique mappings
org_db_genes <- entrez_ids$ENTREZID
mart_genes <- as.character(mart_ids$entrezgene_id)
overlap <- intersect(org_db_genes, mart_genes)

cat("\nMapping comparison:\n")
cat("Genes mapped by both methods:", length(overlap), "\n")
cat("Genes unique to org.Rn.eg.db:", length(setdiff(org_db_genes, mart_genes)), "\n")
cat("Genes unique to biomaRt:", length(setdiff(mart_genes, org_db_genes)), "\n")

# combine the two mapping results
combined_ids <- unique(c(org_db_genes, mart_genes))
combined_ids

#gene_sd <- apply(eset_final, 1, sd)
#variable_genes <- names(gene_sd[gene_sd > 0])
#filtered_matrix <- eset_final[variable_genes, ]

# no gene fulfill the sd > 0 condition
# so filter with top 30% variance

variance_threshold <- 0.2
variance <- apply(eset_final, 1, var)
variable_genes <- names(variance[variance > variance_threshold])
gene_indices <- match(variable_genes, rownames(eset_final))
filtered_matrix <- eset_final[gene_indices, , drop=FALSE]
# Remove rows with NA rownames
filtered_matrix <- filtered_matrix[!is.na(rownames(filtered_matrix)), ]

msigdb_gene_sets <- msigdbr(species = "Rattus norvegicus", category = "H")
gene_sets <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)
# Create parameter object for ssGSEA
params <- ssgseaParam(filtered_matrix, 
                      gene_sets, 
                      minSize = 1,
                      maxSize = Inf)

# Run ssGSEA analysis
ssGSEA_results <- gsva(params)

# Create annotation for samples
annotation_col <- data.frame(
    Group = factor(sapply(strsplit(colnames(ssGSEA_results), "_"), `[`, 1))
)
rownames(annotation_col) <- colnames(ssGSEA_results)

# Generate heatmap
heatmap <- pheatmap(ssGSEA_results,
                    scale = "row",
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    show_colnames = TRUE,
                    show_rownames = TRUE,
                    annotation_col = annotation_col,
                    main = "ssGSEA Pathway Activity Score")

# Save the heatmap
pdf("ssgsea_heatmap.pdf", width = 12, height = 10)
print(heatmap)
dev.off()