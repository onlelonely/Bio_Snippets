# ---------------------------------------------
# Title: MEGENA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/MEGENA.md
# ---------------------------------------------

MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,
                           hub.pval = hub.pval,
                           remove.unsig = TRUE,
                           min.size = modulesize_cutoff,
                           max.size = vcount(g)/2,
                           doPar = doPar,
                           num.cores = n.cores,
                           n.perm = hub.perm,
                           save.output = TRUE)

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,
                                       hub.pvalue = hub.pval,
                                       min.size = modulesize_cutoff,
                                       max.size = vcount(g)/2,
                                       annot.table = NULL,
                                       id.col = NULL,
                                       symbol.col = NULL,
                                       output.sig = TRUE)

module_genes <- summary.output$modules