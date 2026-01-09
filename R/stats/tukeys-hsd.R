# ---------------------------------------------
# Title: Tukey's HSD
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Statistical Tests/Tukey's HSD.md
# ---------------------------------------------

# 假設 'iris' 資料集，我們想比較不同物種 (Species) 的花萼寬度 (Sepal.Width)

# 1. 執行 ANOVA
anova_model <- aov(Sepal.Width ~ Species, data = iris)

# 2. 執行 Tukey's HSD 檢定
tukey_results <- TukeyHSD(anova_model)

# 3. 檢視結果
print(tukey_results)
#           diff        lwr         upr     p adj
# versicolor-setosa    -0.658 -0.8188552 -0.4971448 0.0000000
# virginica-setosa     -0.454 -0.6148552 -0.2931448 0.0000000
# virginica-versicolor  0.204  0.0431448  0.3648552 0.0087802

# 4. 視覺化信賴區間
plot(tukey_results)