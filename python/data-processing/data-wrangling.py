# ---------------------------------------------
# Title: Data Wrangling
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/Data Wrangling.md
# ---------------------------------------------

import pandas as pd
import numpy as np

# 假設 df 是一個 DataFrame

# 清洗：處理缺失值
df['column_name'].fillna(df['column_name'].mean(), inplace=True) # 用平均值填充

# 建構：轉換資料類型
df['date_column'] = pd.to_datetime(df['date_column'])

# 豐富：從現有欄位創建新欄位 (特徵工程)
df['year'] = df['date_column'].dt.year

# 篩選與轉換
df_clean = df[df['value_column'] > 0].copy()
df_clean['log_value'] = np.log(df_clean['value_column'])