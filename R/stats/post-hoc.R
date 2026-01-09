# ---------------------------------------------
# Title: Post-hoc
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Statistical Tests/Post-hoc.md
# ---------------------------------------------

# 假設 data 是一個 data frame，包含數值變數 'value' 和分組變數 'group'

# 1. 執行 ANOVA
anova_result <- aov(value ~ group, data = data)
summary(anova_result)

# 2. 如果 ANOVA 結果顯著，則執行 Tukey's HSD 事後檢定
if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  
  # 視覺化結果
  plot(tukey_result)
}