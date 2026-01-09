# ---------------------------------------------
# Title: GEO analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/GEO analysis.md
# ---------------------------------------------

library(GEOquery)
library(limma)
library(WGCNA)
library(edgeR)

# load series and platform data from GEO

gset <- getGEO("GSE15641", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gsetidx

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11111111111111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXX00000000000XXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")1

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# Multiprobe to single representation
ex <- exprs(gset)
probeToGene <- as.character(unlist(fData(gset)[3]))
collapsedExprData <- collapseRows(ex, rowGroup = probeToGene, rowID = rownames(ex), method = "Average")
ex <- collapsedExprData$datETcollapsed

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("pRCC","Norm"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# Normalization
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 1)] <- NaN
  ex <- log2(ex) }

ex_complete <- ex[complete.cases(ex), ]