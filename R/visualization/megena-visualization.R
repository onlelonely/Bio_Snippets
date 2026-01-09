# ---------------------------------------------
# Title: MEGENA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/MEGENA.md
# ---------------------------------------------

calculate_eigengene <- function(expr_data, genes) {
  # Check if all genes exist in expression data
  valid_genes <- genes[genes %in% rownames(expr_data)]
  
  if (length(valid_genes) == 0) {
    warning("No valid genes found in module")
    return(NULL)
  }
  
  # Use only valid genes for calculation
  module_expr <- expr_data[valid_genes, ]
  
  # If only one gene, return its expression directly
  if (length(valid_genes) == 1) {
    return(as.numeric(module_expr))
  }
  
  # Calculate eigengene using PCA
  pca <- prcomp(t(module_expr), scale. = TRUE)
  return(pca$x[, 1])
}

# Compute eigengenes for each module
module_eigengenes <- lapply(module_genes, function(genes) {
  calculate_eigengene(datExpr, genes)
})

# Convert the list of eigengenes into a dataframe
# Each column is a module, each row is a sample
module_eigengenes_df <- do.call(cbind, module_eigengenes)
colnames(module_eigengenes_df) <- names(module_genes)
rownames(module_eigengenes_df) <- colnames(datExpr)

if ("PredictedResponse" %in% colnames(metadata)) {
  response_vector <- metadata[rownames(module_eigengenes_df), "PredictedResponse"]
  # Run correlations for each module
  cor_results <- apply(module_eigengenes_df, 2, function(ME) {
    test_res <- cor.test(ME, as.numeric(response_vector == "Good"), use = "pairwise.complete.obs", method = "pearson")
    c(correlation = test_res$estimate, p.value = test_res$p.value)
  })
  
  # Format the results
  cor_results_df <- as.data.frame(t(cor_results))
  cor_results_df$Module <- rownames(cor_results_df)
  cor_results_df <- cor_results_df[order(cor_results_df$p.value), ]
} else {
  warning("The specified metadata variable does not exist.")
}

# Create the encoded module names (e.g., M1, M2, ...)
encoded_module_names <- paste0("M", seq_along(cor_results_df$Module))

# Add the encoded module names to the correlation results dataframe
cor_results_df$EncodedModule <- encoded_module_names

# Add the encoded module names to the correlation results dataframe
cor_results_df$EncodedModule <- encoded_module_names

# Calculate the correlation matrix for the module eigengenes
corr_matrix <- cor(module_eigengenes_df, use = "pairwise.complete.obs", method = "pearson")

dist_matrix <- dist(corr_matrix) # Euclidean distance between modules
hc <- hclust(dist_matrix, method = "ward.D2")

# Supermodule 
# Use a height cut to determine the number of super-modules
height_cut <- 20
supermodule_assignment <- cutree(hc, h = height_cut)

# Add super-module assignment to correlation results
cor_results_df$SuperModule <- paste0("SM", supermodule_assignment)

# Combine the data to create the reference dataframe
module_reference <- data.frame(
  Original_Module_Name = cor_results_df$Module,
  Encoded_Module_Name = cor_results_df$EncodedModule,
  SuperModule = cor_results_df$SuperModule,
  stringsAsFactors = FALSE
)

# Save the reference as a CSV file
write.csv(module_reference, "module_reference.csv", row.names = FALSE)

# 1. Aggregate the correlation values by super-module, for example by taking the mean correlation.
supermodule_means <- aggregate(cor_results_df$correlation, 
                               by = list(SuperModule = cor_results_df$SuperModule), 
                               FUN = mean)
colnames(supermodule_means)[2] <- "MeanCorrelation"

# Convert this to a matrix for the heatmap
super_corr_matrix <- as.matrix(supermodule_means$MeanCorrelation)
rownames(super_corr_matrix) <- supermodule_means$SuperModule
colnames(super_corr_matrix) <- "MeanCorrelation"

# Shorten the super-module names if needed
# (They are already fairly short, but you can rename them if desired.)
# For example:
original_sm_names <- rownames(super_corr_matrix)
new_sm_names <- paste0("SuperMod_", seq_along(original_sm_names))
sm_name_ref <- data.frame(
  Original_SuperModule_Name = original_sm_names,
  Short_SuperModule_Name = new_sm_names,
  stringsAsFactors = FALSE
)
write.csv(sm_name_ref, "supermodule_name_reference.csv", row.names = FALSE)
rownames(super_corr_matrix) <- new_sm_names

# Output a heatmap of super-modules
pdf("supermodule_correlation_heatmap.pdf", width = 10, height = 6)
pheatmap(super_corr_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Super-Module Mean Correlation with PredictedResponse",
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Generate a heatmap of the original modules' correlation values
original_corr_matrix <- as.matrix(cor_results_df$correlation)
rownames(original_corr_matrix) <- cor_results_df$Module
colnames(original_corr_matrix) <- "Correlation"

# Output a heatmap of original modules
pdf("original_module_correlation_heatmap.pdf", width = 10, height = 6)
pheatmap(original_corr_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Original Module Correlation with PredictedResponse",
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()