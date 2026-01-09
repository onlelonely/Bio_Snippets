# ---------------------------------------------
# Title: WGCNA
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Bioinformatics Tools/WGCNA.md (4 blocks)
# ---------------------------------------------

# --- Part 1 ---
library(fastDummies)

# Function to generate dummy variables, including interactions
convert_to_dummy_with_interactions <- function(df, column_names) {
  # Create interactions column
  interaction_column <- paste(dfcolumn_names[1], dfcolumn_names[2], sep = "_")
  interaction_column_name <- paste(column_names[1], column_names[2], sep = "_Interaction")
  dfinteraction_column_name <- interaction_column
  
  # Replace NA and blanks with 'No data' for each original column and the interaction column
  for(column_name in c(column_names, interaction_column_name)) {
    dfcolumn_name[is.na(dfcolumn_name)] <- 'No data'
    dfcolumn_name[dfcolumn_name == ''] <- 'No data'
  }
  
  # Generate dummy variables for all selected columns and the interaction column
  df_dummy <- dummy_cols(df, select_columns = c(column_names, interaction_column_name), remove_selected_columns = TRUE)
  
  df_dummy
}

# Apply the function to generate dummies, including for interactions
traits_cat <- convert_to_dummy_with_interactions(traits, c("Drug", "Group"))

rownames(traits_cat) <- rownames(traits)
traits_cat <- traits_cat[,-1]

# --- Part 2 ---
all_traits <- Traits_common
colnames(all_traits) <- c("Age","PSA_BL","PSA_3M","ARB","ASI","Good","Poor","ARB x Good", "ARB x Poor", "ASI x Good","ASI x Poor")

modNames =substring(names(MEs), 3)
geneModuleMembership =as.data.frame(cor(datExpr, MEs, use ="p"));
MMPvalue =as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) =paste("MM", modNames, sep="");
names(MMPvalue) =paste("p.MM", modNames, sep="");
geneTraitSignificance =as.data.frame(cor(datExpr[common_row_names,], all_traits, use ="p"));
GSPvalue =as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) =paste("GS.",names(all_traits), sep="");
names(GSPvalue) =paste("p.GS.",names(all_traits), sep="");

MET = orderMEs(cbind(MEs_common, all_traits))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
pdf("WGCNA-step7-Eigengene-dendrogram.pdf", width = 12, height = 10)
plotEigengeneNetworks(MET, "",excludeGrey = TRUE, marDendro = c(0,4,1,2), marHeatmap = c(6,7,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

moduleColors = mergedColors
module ="midnightblue"
column =match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow =c(1,1));
png("Scatter_MM_GS.png",width = 800,height = 600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab =paste("ModuleMembershipin", module,"module"),
                   ylab ="Gene_significance_for_gleason_score",
                   main =paste("Modulemembershipvs.genesignificance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis= 1.2,col= module)
dev.off()

# --- Part 3 ---
hubs = chooseTopHubInEachModule(datExpr, colorh=moduleColors, power=pow, type=type)
con <- nearestNeighborConnectivity(datExpr, nNeighbors=50, power=pow,
                                   type=type, corFnc = corType)

# --- Part 4 ---
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = pow, verbose = 5)
probes = colnames(datExpr) 
colors_vector <- unlist(strsplit(moduleColors, " "))
unique_colors <- unique(colors_vector)
module = unique_colors
# Select module probes

for (mod in module) {
    inModule = (moduleColors == mod)
    modProbes = probes[inModule]
    modTOM = TOM[inModule, inModule]
    
    dimnames(modTOM) = list(modProbes, modProbes)
    
    cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", mod, ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", mod, ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule]
    )
    # Add any additional code here to process each module
}