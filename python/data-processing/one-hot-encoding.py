# ---------------------------------------------
# Title: One Hot Encoding
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/One Hot Encoding.md
# ---------------------------------------------

import pandas as pd

# 範例 DataFrame
data = {'ID': [1, 2, 3], 'Color': ['Red', 'Green', 'Blue']}
df = pd.DataFrame(data)

# 使用 get_dummies 進行獨熱編碼
df_encoded = pd.get_dummies(df, columns=['Color'], prefix='Color')

print(df_encoded)
#    ID  Color_Blue  Color_Green  Color_Red
# 0   1           0            0          1
# 1   2           0            1          0
# 2   3           1            0          0