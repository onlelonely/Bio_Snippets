# ---------------------------------------------
# Title: TCGA analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/TCGA analysis.md (4 blocks)
# ---------------------------------------------

# --- Part 1 ---
library(survival)
library(survminer)
library(stringr)

#CoxPH                     
cox_results_OS <- data.frame(Gene = character(), Group_HR = numeric(), p = numeric(), stringsAsFactors = FALSE)
cox_results_PFI <- data.frame(Gene = character(), Group_HR = numeric(), p = numeric(), stringsAsFactors = FALSE)
cox_results_OS_5y <- data.frame(Gene = character(), Group_HR = numeric(), p = numeric(), stringsAsFactors = FALSE)
cox_results_PFI_5y <- data.frame(Gene = character(), Group_HR = numeric(), p = numeric(), stringsAsFactors = FALSE)

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Analysize max 5 year survival
merged_survival_normalized$PFI.time_5y<-ifelse(merged_survival_normalized$PFI.time > 5*365, 5*365, merged_survival_normalized$PFI.time) 
merged_survival_normalized$OS.time_5y<-ifelse(merged_survival_normalized$OS.time > 5*365, 5*365, merged_survival_normalized$OS.time) 
merged_survival_normalized$PFI_5y<-ifelse(merged_survival_normalized$PFI.time > 5*365 & merged_survival_normalized$PFI == 1, 0, merged_survival_normalized$PFI)
merged_survival_normalized$OS_5y<-ifelse(merged_survival_normalized$OS.time > 5*365 & merged_survival_normalized$OS == 1, 0, merged_survival_normalized$OS)

# Clean data, remove normal tissue
clean_data_for_survival <- merged_survival_normalized[grep("01[A-B]$", merged_survival_normalized$sample), ]

DEGs <- vector()
# loop through each gene
for (gene in colnames(clean_data_for_survival)[19:91]) {
    # check if the gene is present in the data
    if (sum(colnames(clean_data_for_survival) == gene) > 0) {
        # calculate the median expression level of the gene
        gene_median <- median(clean_data_for_survival[, gene])
        
        # divide the subjects into low and high expression groups based on the median
        clean_data_for_survival$group <- ifelse(clean_data_for_survival[, gene] > gene_median, "High", "Low")
        clean_data_for_survival$group[clean_data_for_survival[, gene] == gene_median] <- ifelse(runif(sum(clean_data_for_survival[, gene] == gene_median)) > 0.5, "Low", "High")
		
		# perform Cox regression for the low and high expression groups
        fit <- coxph(Surv(PFI.time, PFI) ~  relevel(factor(group), ref = "Low"), data = clean_data_for_survival)
        # extract the HRs and p-value and add them to the results data frame
        cox_results_PFI <- rbind(cox_results_PFI, data.frame(Gene = gene, Group_HR = summary(fit)$conf.int[1,1], p = summary(fit)$sctest[3], stringsAsFactors = FALSE))
        
        # Perform the log-rank test
        logrank_test <- survdiff(Surv(PFI.time, PFI) ~ group, data = clean_data_for_survival)
		# Extract the p-value
		p_value_PFI <- pchisq(logrank_test$chisq, length(logrank_test$n) - 1, lower.tail = FALSE)
		if (summary(fit)$conf.int[1,1] > 1.5 & summary(fit)$sctest[3] < 0.05 & p_value_PFI <0.01){
            DEGs <- c(DEGs, gene)
			km <- ggsurvplot(
				  survfit(Surv(PFI.time, PFI) ~ group, data = clean_data_for_survival),
				  data = clean_data_for_survival,
				  pval = TRUE,
				  conf.int = TRUE,
				  ggtheme = theme_bw(),
				  xlab = "Time (days)",
				  ylab = "Probability of survival",
				  legend.title = paste("Gene:", gene),
				  legend.labs = c("High", "Low"),
				  risk.table = TRUE,
				  palette = c("#F8766D", "#00BFC4")
				)
			ggsave(file.path("Output", paste0(gene, "_PFI_km_plot.png")), km, type = "tif", height = 6, width = 6, dpi = 300)
		}		
		
		fit <- coxph(Surv(PFI.time_5y, PFI_5y) ~ relevel(factor(group), ref = "Low"), data = clean_data_for_survival)
        cox_results_PFI_5y <- rbind(cox_results_PFI_5y, data.frame(Gene = gene, Group_HR = summary(fit)$conf.int[1,1], p = summary(fit)$sctest[3], stringsAsFactors = FALSE))
		logrank_test <- survdiff(Surv(PFI.time_5y, PFI_5y) ~ group, data = clean_data_for_survival)
		p_value_PFI_5y <- pchisq(logrank_test$chisq, length(logrank_test$n) - 1, lower.tail = FALSE)	
		if (summary(fit)$conf.int[1,1] > 1.5 & summary(fit)$sctest[3] < 0.05 & p_value_PFI_5y < 0.01){	
			DEGs <- c(DEGs, gene)
            km_5y <- ggsurvplot(
				  survfit(Surv(PFI.time_5y, PFI_5y) ~ group, data = clean_data_for_survival),
				  data = clean_data_for_survival,
				  pval = TRUE,
				  conf.int = TRUE,
				  ggtheme = theme_bw(),
				  xlab = "Time (days)",
				  ylab = "Probability of survival",
				  legend.title = paste("Gene:", gene),
				  legend.labs = c("High", "Low"),
				  risk.table = TRUE,
				  palette = c("#F8766D", "#00BFC4")
				)
			ggsave(file.path("Output", paste0(gene, "_PFI_5y_km_plot.png")), km_5y, type = "tif", height = 6, width = 6, dpi = 300)
		}
	}
}

# loop through each gene
for (gene in colnames(clean_data_for_survival)[19:91]) {
    # check if the gene is present in the data
    if (sum(colnames(clean_data_for_survival) == gene) > 0) {
        # calculate the median expression level of the gene
        gene_median <- median(clean_data_for_survival[, gene])
        
        # divide the subjects into low and high expression groups based on the median
        clean_data_for_survival$group <- ifelse(clean_data_for_survival[, gene] > gene_median, "High", "Low")
        clean_data_for_survival$group[clean_data_for_survival[, gene] == gene_median] <- ifelse(runif(sum(clean_data_for_survival[, gene] == gene_median)) > 0.5, "Low", "High")
		
		# perform Cox regression for the low and high expression groups
        fit <- coxph(Surv(OS.time, OS) ~  relevel(factor(group), ref = "Low"), data = clean_data_for_survival)
        # extract the HRs and p-value and add them to the results data frame
        cox_results_OS <- rbind(cox_results_OS, data.frame(Gene = gene, Group_HR = summary(fit)$conf.int[1,1], p = summary(fit)$sctest[3], stringsAsFactors = FALSE))
        logrank_test <- survdiff(Surv(OS.time, OS) ~ group, data = clean_data_for_survival)
		p_value_OS <- pchisq(logrank_test$chisq, length(logrank_test$n) - 1, lower.tail = FALSE)	
		if (summary(fit)$conf.int[1,1] > 1.5 & summary(fit)$sctest[3] < 0.05 & p_value_OS < 0.01){
			DEGs <- c(DEGs, gene)
            km <- ggsurvplot(
				  survfit(Surv(OS.time, OS) ~ group, data = clean_data_for_survival),
				  data = clean_data_for_survival,
				  pval = TRUE,
				  conf.int = TRUE,
				  ggtheme = theme_bw(),
				  xlab = "Time (days)",
				  ylab = "Probability of survival",
				  legend.title = paste("Gene:", gene),
				  legend.labs = c("High", "Low"),
				  risk.table = TRUE,
				  palette = c("#F8766D", "#00BFC4")
				)
			ggsave(file.path("Output", paste0(gene, "_OS_km_plot.png")), km, type = "tif", height = 6, width = 6, dpi = 300)
		}		
		
		fit <- coxph(Surv(OS.time_5y, OS_5y) ~ relevel(factor(group), ref = "Low"), data = clean_data_for_survival)
        cox_results_OS_5y <- rbind(cox_results_OS_5y, data.frame(Gene = gene, Group_HR = summary(fit)$conf.int[1,1], p = summary(fit)$sctest[3], stringsAsFactors = FALSE))
		logrank_test <- survdiff(Surv(OS.time_5y, OS_5y) ~ group, data = clean_data_for_survival)
		p_value_OS_5y <- pchisq(logrank_test$chisq, length(logrank_test$n) - 1, lower.tail = FALSE)
		if (summary(fit)$conf.int[1,1] > 1.5 & summary(fit)$sctest[3] < 0.05 & p_value_OS_5y <0.01){	
			DEGs <- c(DEGs, gene)
            km_5y <- ggsurvplot(
				  survfit(Surv(OS.time_5y, OS_5y) ~ group, data = clean_data_for_survival),
				  data = clean_data_for_survival,
				  pval = TRUE,
				  conf.int = TRUE,
				  ggtheme = theme_bw(),
				  xlab = "Time (days)",
				  ylab = "Probability of survival",
				  legend.title = paste("Gene:", gene),
				  legend.labs = c("High", "Low"),
				  risk.table = TRUE,
				  palette = c("#F8766D", "#00BFC4")
				)
			ggsave(file.path("Output", paste0(gene, "_OS_5y_km_plot.png")), km_5y, type = "tif", height = 6, width = 6, dpi = 300)
		}
	}
}

matched_genes <- results[results$ensembl_gene_id %in% DEGs,]

# --- Part 2 ---
library(ggplot2)
library(maxstat)

expression_string <- paste("sign(", paste0("(", final_gene_signature$Coefficient, ")*(lasso_data$", final_gene_signature$Gene, ")", collapse = " + "), ") * log(1 + abs(", paste0("(", final_gene_signature$Coefficient, ")*(lasso_data$", final_gene_signature$Gene, ")", collapse = " + "), "))")
lasso_data$risk <- eval(parse(text = expression_string))
lasso_data$risk_score_normalized <- lasso_data$risk - mean(lasso_data$risk)
lasso_data$risk_score_normalized <- lasso_data$risk_score_normalized/sd(lasso_data$risk_score_normalized)

ggplot(lasso_data, aes(x = risk_score_normalized, y = PFI.time/365, color = factor(PFI))) +
    geom_point() +
    labs(x = "Normalized Risk Score", y = "Survival Time", color = "Survival Status")+
    xlim(c(-2,2))+
	scale_color_manual(values = c("1" = "#004B97", "0" = "#D9B300"), labels = c("Progressed", "Not progressed")
    )+ geom_vline(xintercept = -0.54, linetype=4)

PFI_maxstat<- maxstat.test(Surv(as.numeric(lasso_data$PFI.time), lasso_data$PFI) ~ risk_score_normalized, data=lasso_data, smethod="LogRank", pmethod="Lau94")
OS_maxstat<- maxstat.test(Surv(as.numeric(lasso_data$OS.time), lasso_data$OS) ~ risk_score_normalized, data=lasso_data, smethod="LogRank", pmethod="Lau94")

lasso_data$group <- ifelse((lasso_data$risk_score_normalized) <= PFI_maxstat$estimate, "Low", "High")
fit <- survfit(Surv(PFI.time_5y, as.numeric(PFI_5y)) ~ group, data = lasso_data)
p<- ggsurvplot(fit, data = lasso_data,
               palette = c("#F8766D", "#00BFC4"),
               legend.labs = c("High", "Low"),
               legend.title = "PFI 5y",
               risk.table = TRUE,
               pval = TRUE,
               ylim = c(0.3, 1)
               )
ggsave(filename = "riskScore_PFI.png", plot = p, width = 6, height = 6, dpi = 300)

res.logrank <- survdiff(Surv(PFI.time_5y, as.numeric(PFI_5y)) ~ group, data = lasso_data)
pval <- pchisq(res.logrank$chisq, length(res.logrank$n) - 1, lower.tail = FALSE)

# 從merge_data分normal tissue or primary tumor
genes <- c("ENSG00000019549.13","ENSG00000062038.14","ENSG00000092621.13","ENSG00000113721.14","ENSG00000115457.10","ENSG00000145362.20","ENSG00000156298.13","ENSG00000182718.18")
df1_subset <- merged_survival[, 1:27]
df2_subset <- merged_survival[, genes]
prediction_data <- cbind(df1_subset, df2_subset)
prediction_data <- prediction_data %>% 
  mutate_at(vars(28:35), as.numeric)
prediction_data$risk <- sign(((-0.2975)*(prediction_data$ENSG00000019549.13) +
                   (0.0626)*(prediction_data$ENSG00000062038.14) +
                   (-0.0207)*(prediction_data$ENSG00000092621.13) +
                   (-0.0145)*(prediction_data$ENSG00000113721.14) +
                   (-0.0236)*(prediction_data$ENSG00000115457.10) +
                   (-0.0051)*(prediction_data$ENSG00000145362.20) +
                   (0.0831)*(prediction_data$ENSG00000156298.13) +
                   (-0.0015)*(prediction_data$ENSG00000182718.18))) * 
               log(1 + abs(((-0.2975)*(prediction_data$ENSG00000019549.13) +
                   (0.0626)*(prediction_data$ENSG00000062038.14) +
                   (-0.0207)*(prediction_data$ENSG00000092621.13) +
                   (-0.0145)*(prediction_data$ENSG00000113721.14) +
                   (-0.0236)*(prediction_data$ENSG00000115457.10) +
                   (-0.0051)*(prediction_data$ENSG00000145362.20) +
                   (0.0831)*(prediction_data$ENSG00000156298.13) +
                   (-0.0015)*(prediction_data$ENSG00000182718.18)))) 
prediction_data$risk_score_normalized <- prediction_data$risk - mean(prediction_data$risk)
prediction_data$risk_score_normalized <- prediction_data$risk_score_normalized/sd(prediction_data$risk_score_normalized)

prediction_data_tumor <-subset(prediction_data, grepl("01[A-B]$", sample))

PFI_maxstat<- maxstat.test(Surv(as.numeric(prediction_data_tumor$PFI.time), prediction_data_tumor$PFI) ~ risk_score_normalized, data=prediction_data_tumor, smethod="LogRank", pmethod="Lau94")

High_risk <- subset(prediction_data_tumor, risk_score_normalized <= PFI_maxstat$estimate)
Low_risk <- subset(prediction_data_tumor, risk_score_normalized > PFI_maxstat$estimate)

# --- Part 3 ---
library(forestploter)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

# Multicovarite CoxPH
fit_cox <- coxph(Surv(PFI.time,PFI) ~ age_at_initial_pathologic_diagnosis
                 + relevel(factor(TNM_cat), ref = "1") 
                 + relevel(factor(group), ref="Low") 
                 + as.numeric(stage_event_psa) 
                 + relevel(factor(gleason_stage_cat), ref = "<8"), data = lasso_data)
summary(fit_cox)
fit_summary <- summary(fit_cox)

# Plot the multicovariate CoxPH
plot_cox <- as.data.frame(fit_summary$coefficients)
rownames(plot_cox) <- c('Age at diagnosis', 'TNM stage', 'gene signature', 'PSA', 'Gleason stage')
plot_cox$Factor <- rownames(plot_cox)
plot_cox$HR <- plot_cox$`exp(coef)`
plot_cox$upper <-plot_cox$`exp(coef)` + 1.96*plot_cox$`se(coef)`
plot_cox$lower <-plot_cox$`exp(coef)` - 1.96*plot_cox$`se(coef)`
plot_cox <- plot_cox[order(-plot_cox$HR),]
plot_cox$combined <- paste(round(plot_cox$HR, 2), " ", "(", round(plot_cox$lower, 2), " - ", round(plot_cox$upper, 2), ")", sep="")
plot_cox$t <- "                       "
colnames(plot_cox)[11] <-""
colnames(plot_cox)[10] <-"HR (95% CI)"

p<- forest(plot_cox[,c(6,10,11)],
           est = plot_cox$HR,
           lower = plot_cox$lower,
           upper = plot_cox$upper,
           sizes = plot_cox$`se(coef)`,
           ref_line = 1,
           ci_column = 3,
           xlim = c(0, 5),
           ticks_at = c(0, 1, 2, 3),
	       arrow_lab = c("PFI longer", "PFI shorter"),
)
plot(p)

# --- Part 4 ---
library(survival)
library(survminer)
library(stringr)
library(dplyr)
library(ggplot2)
library(timeROC)
library(pROC)

genes <- c("ENSG00000001631","ENSG00000145996","ENSG00000127995","ENSG00000118564")
# Initialize empty list to store coefficient values
coef_values <- list()
for(i in 1:length(gene_signature_coef)){
  # Extract rownames (gene names)
  gene_names <- rownames(gene_signature_coefi)
  # Find the indices of the genes of interest
  gene_indices <- which(gene_names %in% genes)
  # Extract the coefficients for these indices
  gene_coefs <- gene_signature_coefi[gene_indices, ]
  # If there is any 0 in gene_coefs, skip this iteration
  if(any(gene_coefs == 0)) {
    next
  }
  # Append the coefficients to our list
  coef_valuesi <- gene_coefs
}
# Convert list to matrix
coef_matrix <- do.call(rbind, coef_values)

# Calculate column means
average_coefs <- colMeans(coef_matrix, na.rm = TRUE)

roc_curves <- list()
#fit_cox <- coxph(Surv(PFI.time_5y,PFI_5y) ~ age_at_initial_pathologic_diagnosis + relevel(factor(gender), ref ="FEMALE") + relevel(factor(TNM_cat), ref = "1") + relevel(factor(group), ref="Low") + relevel(factor(platelet_qualitative_result), ref = "Normal"), data = lasso_data)

for(gene_id in genes) {
  # Extract the marker values for the current gene
  marker_values <- -lasso_datagene_id  
  # Compute the ROC curve using timeROC
  roc_curvesgene_id <- timeROC(T = lasso_data$PFI.time /365,
                                   delta = lasso_data$PFI,
                                   marker = marker_values,
                                   cause = 1,
                                   weighting = "marginal",
                                   times = c(1, 3, 5),
                                   iid = TRUE)
}


par(pty="s") 

plot(roc_curves"ENSG00000001631", time=3, col="Red", title="")
plot(roc_curves"ENSG00000145996", time=3, add=TRUE, col="Blue")
plot(roc_curves"ENSG00000127995", time=3, add=TRUE, col="Green")
plot(roc_curves"ENSG00000118564", time=3, add=TRUE, col="Purple")
legend("bottomright",col=c("red","blue","green","purple"), 
       legend = c("KRIT1","CDKAL1","CASD1","FBXL5"), lty = 1, bty = "n", cex = 0.8)