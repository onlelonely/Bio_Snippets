# ---------------------------------------------
# Title: Res_pRCC論文revision
# Description: From: Source/3. Efforts/_Archives/Res_pRCC論文revision.md
# ---------------------------------------------

library(GEOquery)
library(Limma)
library(umap)

gset <- getGEO("GSE7023", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL4866", attr(gset, "names")) else idx <- 1
gset <- gsetidx

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000000000000000000000X000000111X11111111"
sml <- strsplit(gsms, split="")1

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Sample","Control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

# fit linear model
fit <- lmFit(exprs(gset), design)
# define contrasts
contrast.matrix <- makeContrasts(Control-Sample, levels=design)
# compute contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
# apply empirical Bayes statistics
fit2 <- eBayes(fit2)
# create top table of results
tT <- topTable(fit2, number=nrow(fit2), adjust="BH")
# add gene symbols to the results
tT <- cbind(Gene.symbol=fData(gset)[rownames(tT), "Gene.symbol"], tT)