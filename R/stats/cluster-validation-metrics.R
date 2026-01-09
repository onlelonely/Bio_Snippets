# ---------------------------------------------
# Title: Clustering
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/Clustering.md
# ---------------------------------------------

# Function to calculate average silhouette width
library(cluster)
library(fpc)
calc_sil <- function(k) {
  clusters <- cutree(hc, k=k)
  sil <- silhouette(clusters, dist_matrix)
  mean(sil[,3])
}

# Function to calculate Calinski-Harabasz index
calc_ch <- function(k) {
  clusters <- cutree(hc, k=k)
  ch <- calinhara(dist_matrix, clusters)
  return(ch)
}

# Test different k values
k_range <- 2:15
sil_scores <- sapply(k_range, calc_sil)
ch_scores <- sapply(k_range, calc_ch)

# Plot results
par(mfrow=c(2,1))
plot(k_range, sil_scores, type="b", pch=19, 
     xlab="Number of clusters (k)", ylab="Average Silhouette Width",
     main="Silhouette Method")
abline(v=k_range[which.max(sil_scores)], lty=2, col="red")

plot(k_range, ch_scores, type="b", pch=19,
     xlab="Number of clusters (k)", ylab="Calinski-Harabasz Index",
     main="Calinski-Harabasz Method")
abline(v=k_range[which.max(ch_scores)], lty=2, col="red")

# Print optimal k values
cat("Optimal k according to Silhouette method:", k_range[which.max(sil_scores)], "\n")
cat("Optimal k according to Calinski-Harabasz method:", k_range[which.max(ch_scores)], "\n")