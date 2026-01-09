# ---------------------------------------------
# Title: WGCNA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/WGCNA.md
# ---------------------------------------------

library(WGCNA)
library(DESeq2)
library(Biomart)

m.vars = apply(eset_final, 1, var)
lower_bound = quantile(m.vars, probs = 0.55) # Top 5=45%
upper_bound = quantile(m.vars, probs = 0.85) # Top 15%
mrna_df_selected = eset_final[which(m.vars > lower_bound & m.vars <= upper_bound), ]

datExpr=as.data.frame(t(mrna_df_selected))
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
table(clust)
## no outlier

# soft threshold
datExpr=as.data.frame(t(mrna_df_selected))
powers = c(c(1:15))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.70,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

pow=8
datExpr=as.data.frame(t(mrna_df_selected))
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

cor <- WGCNA::cor
type = "signed"
corType = "bicor"
net = blockwiseModules(datExpr, power = pow, maxBlockSize = 12000,
                       TOMType = "signed Nowick", minModuleSize = 25,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       verbose =3)
table(net$colors)
cor<-stats::cor
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms1, mergedColors[net$blockGenes1],
                    groupLabels = c("Module colors", 
                                    "GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
					
moduleColors=mergedColors
					
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

common_row_names <- intersect(rownames(MEs), rownames(traits_cat))
MEs_common <- MEs[common_row_names, ]
Traits_common <- traits_cat[common_row_names, ]

moduleTraitCor = cor(MEs_common, Traits_common, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor) 
pdf(file = 'cormap_nowick.pdf',width = 15,height = 15)
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(Traits_common),
               yLabels = names(MEs_common), 
               ySymbols = names(MEs_common),
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.7, 
               zlim = c(-1,1), 
               cex.lab.y = 0.6,
               cex.lab.x = 0.7,
               yColorWidth = strheight("M")/2,
               main = paste("Module-trait relationships"))
dev.off()[]