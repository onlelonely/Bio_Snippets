# ---------------------------------------------
# Title: Frequency Encoding
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/Frequency Encoding.md
# ---------------------------------------------

import pandas as pd

# 範例 DataFrame
data = {'city': ['Taipei', 'New York', 'Taipei', 'Tokyo', 'New York', 'Taipei']}
df = pd.DataFrame(data)

# 1. 計算頻率並建立映射字典
frequency_map = df['city'].value_counts().to_dict()
# frequency_map -> {'Taipei': 3, 'New York': 2, 'Tokyo': 1}

# 2. 轉換欄位
df['city_freq_encoded'] = df['city'].map(frequency_map)

print(df)