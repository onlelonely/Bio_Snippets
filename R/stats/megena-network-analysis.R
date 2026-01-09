# ---------------------------------------------
# Title: MEGENA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/MEGENA.md
# ---------------------------------------------

# From TC probe selection

library(MEGENA)
library(igraph)
library(DESeq2)
library(Biomart)
library(reshape2)

# Set parameters for MEGENA analysis
density_cutoff = 0.25
modulesize_cutoff = 10
hub_gene_cutoff = 0.1
method = "pearson"
module.pval = 0.05
hub.pval = 0.05
core.perm = 10
hub.perm = 100
n.cores = parallel::detectCores() - 4


# Select probes based on variance
m.vars = apply(eset_final, 1, var)
lower_bound = quantile(m.vars, probs = 0.6) # Top 1-N
upper_bound = quantile(m.vars, probs = 0.95) 
mrna_df_selected = eset_final[which(m.vars > lower_bound & m.vars <= upper_bound), ]

datExpr = as.data.frame(mrna_df_selected)
ijw <- calculate.correlation(datExpr,doPerm = core.perm,
                             output.corTable = FALSE,output.permFDR = FALSE)
# Further parameters
edgelist = ijw[,1:3]
max.skipEdges = NULL
maxENum = NULL
doPar = TRUE
keep.track = TRUE

el <- calculate.PFN(edgelist, doPar = doPar, num.cores = n.cores)
g <- graph.data.frame(el, directed = FALSE)