# ---------------------------------------------
# Title: TCGA analysis
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Databases/TCGA analysis.md (3 blocks)
# ---------------------------------------------

# --- Part 1 ---
query_clin <- GDCquery(project = "TCGA-BRCA", 
				data.category = "Clinical", 
				data.format = "BCR XML", 
				data.type = "Clinical Supplement" )

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type ="STAR - Counts",
                  sample.type = c("Primary Tumor","Solid Tissue Normal")
                 )

query_mutation <- GDCquery(project = "TCGA-BRCA", 
					data.category = "Simple Nucleotide Variation", 
					data.type = "Masked Somatic Mutation", 
					workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking" 
					)

# --- Part 2 ---
GDCdownload(
  query,
  method = "api",
  files.per.chunk = 100,
  directory = "./Data"
)

# --- Part 3 ---
#Filtered with gene list from GO/KEGG
genes <- genes <- c("FTO","CDKAL1","FBXL5","NDUFB4","ITGBL1","SNRPB",
                    "KRIT1","NAA35","NAA15","SMC3","CASD1","USP16","SPAG9",
                    "EGLN1","CD46","SPOPL")
rownames(sig.genes) <- gsub("\\..*", "", rownames(sig.genes))
filtered_mapping <- results[results$hgnc_symbol %in% genes, ]
ensembl_ids <- filtered_mapping$ensembl_gene_id
filtered_sig_genes <- sig.genes[rownames(sig.genes) %in% ensembl_ids, ]
mrna_df_survival <- mrna_df_filtered[cleaned_rownames %in% rownames(filtered_sig_genes), ]

# Create Survival dataset PRAD
sampleMap <- data.frame(sample = colnames(mrna_df_survival), bcr_patient_barcode = substring(colnames(mrna_df_survival), 1, 12))
survival_data_1 <- subset(clinical_extend, select = c("bcr_patient_barcode",
                                                      "age_at_initial_pathologic_diagnosis", 
                                                      "vital_status", 
                                                      "OS", "OS.time",
                                                      "DSS","DSS.time",
                                                      "DFI", "DFI.time", 
                                                      "PFI", "PFI.time"))
survival_data_2 <- subset(clinical_data, select = c("bcr_patient_barcode","stage_event_psa", "stage_event_gleason_grading", "stage_event_tnm_categories"))
survival_data <- merge(survival_data_1, survival_data_2, by = "bcr_patient_barcode")
merged_data <- distinct(merge(sampleMap, survival_data, by = "bcr_patient_barcode"))

# Create Survival dataset BLCA
sampleMap <- data.frame(sample = colnames(mrna_df_survival), bcr_patient_barcode = substring(colnames(mrna_df_survival), 1, 12))
colnames(clinical_extend)[2] <- "bcr_patient_barcode"
survival_data_1 <- subset(clinical_extend, select = c("bcr_patient_barcode", "OS", "OS.time","DSS","DSS.time","DFI", "DFI.time", "PFI", "PFI.time"))
survival_data_2 <- subset(clinical_data, select = c("bcr_patient_barcode","gender","diagnosis_subtype", "neoplasm_histologic_grade","age_at_initial_pathologic_diagnosis", "lymphovascular_invasion_present", "stage_event_pathologic_stage", "stage_event_tnm_categories"))
survival_data <- merge(survival_data_1, survival_data_2, by = "bcr_patient_barcode")
merged_data <- distinct(merge(sampleMap, survival_data, by = "bcr_patient_barcode"))

#convert gleason score (PRAD only)
merged_data$stage_event_gleason_grading <- ifelse(substr(merged_data$stage_event_gleason_grading, 1, 1) %in% c("6", "7", "8", "9"), 
                                                 substr(merged_data$stage_event_gleason_grading, 1, 1), 
                                                 substr(merged_data$stage_event_gleason_grading, 1, 2))

## gleason grouping
first_digit <- as.numeric(substr(merged_data$stage_event_gleason_grading, 1, 1))

# Apply conditions to categorize gleason_stage_cat
merged_data$gleason_stage_cat <- ifelse(
  first_digit == 1, ">=8",
  ifelse(first_digit < 8, "<8", ">=8")
)

### TNM grouping 
# Create an empty vector to store the group assignments
TNM_groups <- numeric(length(merged_data$stage_event_tnm_categories))
# Loop through each subject's TMN data and assign them to a group based on the pattern
for (i in 1:length(merged_data$stage_event_tnm_categories)) {
  if (str_detect(as.character(merged_data$stage_event_tnm_categories[i]), "^(?:(?!.*T[3-9]).*)?(?=.*[T][12]).*(?=.*[M]?0).*[N]?0.*$")) {
	if (str_detect(as.character(merged_data$stage_event_tnm_categories[i]), "^(?=.*[T][3-9]).*$")){
		TNM_groups[i] <- 2
	}
	else {
		TNM_groups[i] <- 1
	}
  } else {
    TNM_groups[i] <- 2
  }
}

# Add the group assignments as a new column to the data frame
merged_data$TNM_cat <- ifelse(TNM_groups == 1, "1", "2")

###Get Normalized Gene expression by vooma & merge
voom_normalized_data <- v$E
rownames(voom_normalized_data) <- c(gsub("\\..*", "", rownames(voom_normalized_data)))
filtered_voom <- voom_normalized_data[rownames(voom_normalized_data) %in% filtered_mapping$ensembl_gene_id,]
normalized_gene_expression_df <- t(as.data.frame(filtered_voom))
row_names <- row.names(normalized_gene_expression_df)
normalized_gene_expression_df <- cbind(sample = row_names, normalized_gene_expression_df)
merged_survival_normalized <- merge(merged_data, normalized_gene_expression_df, by = "sample")