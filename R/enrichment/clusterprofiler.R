# ---------------------------------------------
# Title: ClusterProfiler
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/ClusterProfiler.md
# ---------------------------------------------

# Run GO enrichment analysis
go_enrich <- enrichGO(gene = combined_ids,
                      OrgDb = org.Rn.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Run KEGG pathway analysis
kegg_enrich <- enrichKEGG(gene = combined_ids,
                         organism = "rno",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)


# Create dotplot visualizations
dotplot(go_enrich, showCategory=15, title="GO Enrichment")
dotplot(kegg_enrich, showCategory=15, title="KEGG Pathway Enrichment")

# Create enrichment map for GO terms
# emapplot(pairwise_termsim(go_enrich), showCategory = 30)

# save the plots
pdf("go_enrichment_dotplot.pdf", width = 10, height = 8)
dotplot(go_enrich, showCategory=15, title="GO Enrichment")
dev.off()

pdf("kegg_enrichment_dotplot.pdf", width = 10, height = 8)
dotplot(kegg_enrich, showCategory=15, title="KEGG Pathway Enrichment")
dev.off()