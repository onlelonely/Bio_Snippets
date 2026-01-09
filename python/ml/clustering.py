# ---------------------------------------------
# Title: Clustering
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/Clustering.md (2 blocks)
# ---------------------------------------------

# --- Part 1 ---
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_moons
import matplotlib.pyplot as plt

# Generate sample data
X, y = make_moons(n_samples=200, noise=0.05, random_state=0)

# Fit DBSCAN
dbscan = DBSCAN(eps=0.3, min_samples=5)
clusters = dbscan.fit_predict(X)

# Plot results
plt.scatter(X[:, 0], X[:, 1], c=clusters, cmap="viridis", marker="o")
plt.title("DBSCAN Clustering")
plt.xlabel("Feature 1")
plt.ylabel("Feature 2")
plt.show()

# --- Part 2 ---
from sklearn.mixture import GaussianMixture
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt

# Generate sample data
X, y = make_blobs(n_samples=300, centers=4, cluster_std=0.60, random_state=0)

# Fit GMM
gmm = GaussianMixture(n_components=4, random_state=0)
clusters = gmm.fit_predict(X)

# Plot results
plt.scatter(X[:, 0], X[:, 1], c=clusters, cmap="viridis", marker="o")
plt.title("Gaussian Mixture Model Clustering")
plt.xlabel("Feature 1")
plt.ylabel("Feature 2")
plt.show()