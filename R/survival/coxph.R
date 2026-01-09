# ---------------------------------------------
# Title: CoxPH
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Statistical Tests/CoxPH.md
# ---------------------------------------------

surv_obj <- Surv(time = data$time, event = data$status)
fit.coxph <- coxph(surv_obj ~ age + sex + smoke, data = data)