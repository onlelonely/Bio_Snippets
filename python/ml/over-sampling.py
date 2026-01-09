# ---------------------------------------------
# Title: Over-sampling
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/Over-sampling.md
# ---------------------------------------------

from collections import Counter
from sklearn.datasets import make_classification
from imblearn.over_sampling import SMOTE

# 產生一個不平衡的資料集
X, y = make_classification(n_classes=2, class_sep=2,
                           weights=[0.1, 0.9], n_informative=3, n_redundant=1, flip_y=0,
                           n_features=20, n_clusters_per_class=1, n_samples=1000, random_state=10)

print(f'Original dataset shape %s' % Counter(y))
# Original dataset shape Counter({1: 900, 0: 100})

# 使用 SMOTE 進行上採樣
sm = SMOTE(random_state=42)
X_res, y_res = sm.fit_resample(X, y)

print(f'Resampled dataset shape %s' % Counter(y_res))
# Resampled dataset shape Counter({1: 900, 0: 900})