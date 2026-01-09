# ---------------------------------------------
# Title: Clustering
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/Clustering.md
# ---------------------------------------------

library(ConsensusClusterPlus)
library(ggplot2)

# Filter top 70% variance
temp <- t(eset_final)
data_unique <- temp[!duplicated(temp), ] 
variances <- apply(t(data_unique), 1, var) # Calculate variance for each gene
threshold <- quantile(variances, 0.3) # Find the 730th percentile
high_variance_genes <- names(variances[variances > threshold])
filtered_data <- t(data_unique)[high_variance_genes, ]

# Consensus Clustering on high-variance data
filtered_data <- sweep(filtered_data, 1, apply(filtered_data, 1, median, na.rm=T))
title <- "High_variance_baseline"
results <- ConsensusClusterPlus(filtered_data, maxK=5, reps=100, pItem=0.8, pFeature=0.8, 
                                title=title, clusterAlg="hc", distance="pearson", plot="png")

# Calculate and plot Item Consensus
icl <- calcICL(results, title, plot="png")

# Extract consensus classes
consensus_classes <- results4$consensusClass
write.csv(consensus_classes, "HighVarCcluster.csv")

# Perform PCA

pca_res <- prcomp(t(eset_final), scale. = TRUE)
# Extract the PCA scores
pca_scores <- pca_res$x
# Combine PCA scores with cluster information
pca_data <- data.frame(pca_scores, cluster = consensus_classes)

# Plot PCA with clusters
ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(cluster))) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data", x = "PC1", y = "PC2", color = "Cluster")

library(Rtsne)
tsne_res <- Rtsne(data_unique, dims = 2, Perplexity = 20, verbose = TRUE, max_iter = 3000)
# Extract the t-SNE scores
tsne_scores <- tsne_res$Y
# Combine t-SNE scores with cluster information
tsne_data <- data.frame(tsne_scores, cluster = consensus_classes)

ggplot(tsne_data, aes(x = X1, y = X2, color = factor(cluster))) +
  geom_point() +
  theme_minimal() +
  labs(title = "t-SNE of Gene Expression Data", x = "t-SNE 1", y = "t-SNE 2", color = "Cluster")