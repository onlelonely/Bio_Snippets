# ---------------------------------------------
# Title: BiomaRt
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/BiomaRt.md
# ---------------------------------------------

library(biomaRt)

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# ensembl = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", mirror= "asia")

results = getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'), 
  filters = 'hgnc_symbol',
  values = entrez_ids,
  mart = ensembl
)