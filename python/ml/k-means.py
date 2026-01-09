# ---------------------------------------------
# Title: K-means
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/K-means.md
# ---------------------------------------------

from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt

# 產生範例數據
X, _ = make_blobs(n_samples=300, centers=4, cluster_std=0.60, random_state=0)

# 初始化並擬合K-means模型
# n_clusters 即為 K 值
kmeans = KMeans(n_clusters=4, random_state=0, n_init=10) # n_init=10 避免局部最優
kmeans.fit(X)

# 獲取預測的群集標籤和質心
y_kmeans = kmeans.predict(X)
centers = kmeans.cluster_centers_

# 視覺化結果
plt.scatter(X[:, 0], X[:, 1], c=y_kmeans, s=50, cmap='viridis')
plt.scatter(centers[:, 0], centers[:, 1], c='red', s=200, alpha=0.75, marker='X')
plt.title("K-Means Clustering")
plt.show()