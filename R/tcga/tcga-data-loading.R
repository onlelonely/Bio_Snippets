# ---------------------------------------------
# Title: TCGA analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/TCGA analysis.md
# ---------------------------------------------

## 以PRAD harmonized dataset 為例
library(TCGAbiolinks)
library(dplyr)
library(stringr) 
library(ggplot2)
library(SummarizedExperiment)

### Load functions for identifying upregulated and downregulated genes
get_upregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange >= 1)],
                   rownames(df)[which(df$pvalue <= 0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key), ])
  return(results)
}

get_downregulated <- function(df) {
  key <- intersect(rownames(df)[which(df$log2FoldChange <= -1)],
                   rownames(df)[which(df$pvalue <= 0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key), ])
  return(results)
}

### Load data
mrna_df <- GDCprepare(query, directory = "./data")
clinical_data <- GDCprepare_clinic(query_clin, clinical.info = "patient", directory = "./data")
mutation_data <- GDCprepare(query_mutation, directory = "./data")                    
clinical_extend <- read.csv("./data/PRAD_clinical_extend.csv") #從別人的paper來的survival資料

### Process mutation data to identify TP53 mutations 
tp53_mutations <- mutation_data[mutation_data$Hugo_Symbol == "TP53", ] tp53_mutated_samples <- unique(tp53_mutations$Tumor_Sample_Barcode)

### Function to extract patient ID from barcode 
extract_patient_id <- function(barcode) { paste(unlist(strsplit(barcode, "-"))[1:3], collapse = "-") } 
# Get patient IDs with TP53 mutations 
tp53_mutated_patients <- unique(sapply(tp53_mutated_samples, extract_patient_id))

### Modify column names, set design to specific gene mutation (here is TP53)
mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_df <- assay(mrna_df)

delim_fn = function(x, n, i) {
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(
      X, "-"
    ))[(n + 1):(i)], collapse = "-")))
}

# Extract patient IDs from expression data
colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 3) # Extract first 3 parts for patient ID
                    
mrna_meta <- as.data.frame(mrna_meta)
mrna_df_plain <- as.data.frame(mrna_df)

### Create TP53 mutation status for each sample
mrna_meta$patient_id <- colnames(mrna_df)
mrna_meta$TP53_status <- ifelse(mrna_meta$patient_id %in% tp53_mutated_patients, 
                               "TP53_Mutated", 
                               "TP53_Wild_Type")

### Clean up metadata
mrna_meta <- mrna_meta[, c("patient_id", "TP53_status")]
colnames(mrna_meta) <- c("cases", "Condition")
mrna_meta$Condition <- as.factor(mrna_meta$Condition)

### Filter to keep only samples with clear TP53 status
# Remove any samples that couldn't be classified
mrna_meta <- mrna_meta[!is.na(mrna_meta$Condition), ]
mrna_df_plain <- mrna_df_plain[, mrna_meta$cases]

print(paste("Total samples for analysis:", ncol(mrna_df_plain)))
print(paste("TP53 mutated:", sum(mrna_meta$Condition == "TP53_Mutated")))
print(paste("TP53 wild-type:", sum(mrna_meta$Condition == "TP53_Wild_Type")))

### ------------------------------------------- ###

### modify column names, set design to tumor vs normal
mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_df <- assay(mrna_df)
delim_fn = function(x, n, i) {
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(
      X, "-"
    ))[(n + 1):(i)], collapse = "-")))
}
colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)
                    
mrna_meta <- as.data.frame(mrna_meta)
mrna_df_plain <- as.data.frame(mrna_df)
mrna_meta[, 2] <- as.character(mrna_meta[, 2])
mrna_meta[, 2] <-
  gsub("Primary solid Tumor", "Tumor", mrna_meta[, 2])
mrna_meta[, 2] <-
  gsub("Solid Tissue Normal", "Normal", mrna_meta[, 2])
mrna_meta[, 2] <- as.factor(mrna_meta[, 2])
colnames(mrna_meta) <- c("cases", "Condition")