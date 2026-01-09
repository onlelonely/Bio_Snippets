# ---------------------------------------------
# Title: Hyperparameter
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/Hyperparameter.md
# ---------------------------------------------

from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification

# 準備資料
X, y = make_classification(n_samples=1000, n_features=20, random_state=42)

# 定義模型
rf = RandomForestClassifier()

# 定義超參數搜索空間
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, None],
    'min_samples_leaf': [1, 2, 4]
}

# 設置網格搜索
grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=3, n_jobs=-1, verbose=2)

# 執行搜索
grid_search.fit(X, y)

# 輸出最佳參數
print(f"Best parameters found: {grid_search.best_params_}")