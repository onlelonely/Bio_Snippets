# ---------------------------------------------
# Title: KNN
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Machine Learning/KNN.md
# ---------------------------------------------

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris

# 載入數據
iris = load_iris()
X, y = iris.data, iris.target

# 分割數據集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 特徵縮放
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 初始化並擬合KNN模型 (選擇K=3)
knn = KNeighborsClassifier(n_neighbors=3)
knn.fit(X_train_scaled, y_train)

# 進行預測與評估
accuracy = knn.score(X_test_scaled, y_test)
print(f"Accuracy: {accuracy:.4f}")