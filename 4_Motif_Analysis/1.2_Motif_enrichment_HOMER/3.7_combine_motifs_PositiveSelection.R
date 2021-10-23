# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ComplexHeatmap) # For Heatmap function
library(circlize) # for generating colors using colorremp2
library(gplots) # for heatmap.2 which has better row names
library (dplyr)

# Function to read motif file and get selected columns ---------------------------------------------------------------------------------------
# Input: motif file name and a prefic for column name
# Output: motif data frame with selected columns
ReadMotifFile = function(motiffile, col.prefix) {
  
  df = read.table(file=motiffile, head = T, sep = ",")
  #df = df[, c(1,3,5,7,9, 10:14)]
  #colnames(df) = c(paste0(col.prefix, ".motif"), paste0(col.prefix, ".pval"), paste0(col.prefix, ".fdr"), paste0(col.prefix, ".fg"), paste0(col.prefix, ".bg"), paste0(col.prefix, ".FC"), paste0(col.prefix, ".logFC"), paste0(col.prefix, ".Escore"), paste0(col.prefix, ".clusterId"), "finalGene")
  df = df[, c(3,5,7,9, 10:14)]
  colnames(df) = c(paste0(col.prefix, ".pval"), paste0(col.prefix, ".fdr"), paste0(col.prefix, ".fg"), paste0(col.prefix, ".bg"), paste0(col.prefix, ".FC"), paste0(col.prefix, ".logFC"), paste0(col.prefix, ".Escore"), paste0(col.prefix, ".clusterId"), "finalGene")
  
  df =  df %>% 
    group_by(finalGene) %>% 
    summarise(across(everything(), list(min)))
  
  return(df)
}

# Read motif enrichment file ------------------------------------------------------------------------------------------------------------------
# Nfur all UP PEAKS
nfur.all = ReadMotifFile("Motif_enrichment_results/Nfur_Master_Up/Nfur_All/knownResults_processed.csv", "nfur.all")
head(nfur.all)

# Nfur positive selection peaks
posSel = ReadMotifFile("Motif_enrichment_results/Positive_selection/knownResults_processed.csv", "nfur.pos.sel")
head(posSel)

merged <- Reduce(function(x, y) merge(x, y, by="finalGene", all = TRUE), list(nfur.all, posSel))
head(merged)

merged  = as.data.frame(merged)
rownames(merged) = merged$finalGene

# Select the motifs significant in at least one conditions
toPlot = merged[(merged$nfur.all.fdr_1 <= 0.1 |
                 merged$nfur.pos.sel.fdr_1 <= 0.1
                          ),]
head(toPlot)
colnames(head(toPlot))

# Now get 3 sets of genes nfur, killifish and conserved and get a file sorted based on ages ----------------------------------------
# Only significant in all motif enrichment
allOnly = toPlot[(toPlot$nfur.all.fdr_1 <= 0.1 &
                  toPlot$nfur.pos.sel.fdr_1 > 0.1
                         ),]
# Only significant in positive selection motif enrichment
positiveOnly = toPlot[(toPlot$nfur.all.fdr_1 > 0.1 &
                       toPlot$nfur.pos.sel.fdr_1 <= 0.1
),]

# Significant in both
both = toPlot[(toPlot$nfur.all.fdr_1 <= 0.1 &
                 toPlot$nfur.pos.sel.fdr_1 <= 0.1
),]

# Below condition must be true
sum (length(allOnly[,1]), length(positiveOnly[,1]), length(both[,1])) == length(toPlot[,1])

# Merge and save in anew file
toPlotNew = rbind(both, positiveOnly, allOnly)

write.csv(toPlotNew, file = "Output_files_merged/positiveSelectionProcessed.csv")
write.csv(toPlot, file = "Output_files_merged/positiveSelectionFinal.csv")

# Mark only the top motifs to be labelled on the plot. These are selected significant labels that we want to highlight
genesLabels = toPlotNew$finalGene
index = match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARG(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels)
match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARA(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels) # Nothing should be NA
custom_label_genes = rowAnnotation(foo = anno_mark(at = index,
                                                   side = "left", # To put labels on left or riht side
                                   labels = genesLabels[index]))

# Print file with labels
pdf(paste0("Output_plots/", "positive_selection", "_selectedLabs.pdf"), width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11)]), 
        col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "grey80")), 
        cluster_rows = F, 
        cluster_columns = F,
        show_row_names = F,
        left_annotation = custom_label_genes # To put labels on left or riht side
)
dev.off()

# Plot another heatmap with all labels 
pdf(paste0("Output_plots/", "positive_selection", "_AllLabs.pdf"), width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11)]),
        col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "grey80")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 7),
        row_names_max_width = unit(30, "cm"))
dev.off()

