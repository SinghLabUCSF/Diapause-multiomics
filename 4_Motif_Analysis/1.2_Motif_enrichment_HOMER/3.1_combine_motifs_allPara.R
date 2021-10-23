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

# Run this for 1bp, 25% or 50% overlap *************************************************************
rm (list = setdiff(ls(), "ReadMotifFile"))
overlap = 'Nfur_UP_01' # 'Nfur_UP_01' 'Nfur_UP_25' 'Nfur_UP_50'   *** Peaks with 1bp, 25% or 50% overlap
peaks = 'Either'                                              #   *** genomic regions with or without peaks

preFix = paste("NeoF", overlap, sep = "_")                    #   *** File name prefix - NeoF combined dataset
preFix

# Read motif enrichment files ------------------------------------------------------------------------------------------------------------------
# Nfur all UP PEAKS
nfur.all = ReadMotifFile("Motif_enrichment_results/Nfur_Master_Up/Nfur_Neo_Combine/knownResults_processed.csv", "nfur.all")
head(nfur.all)

aaul.all = ReadMotifFile(paste0("Motif_enrichment_results/Nfur_Master_Up/", overlap, "/" , peaks, "/Aaus_Neo_Combine/knownResults_processed.csv"), "aaul.all")
head(aaul.all)

ast.all = ReadMotifFile(paste0("Motif_enrichment_results/Nfur_Master_Up/", overlap, "/" , peaks, "/Astr_Neo_Combine/knownResults_processed.csv"), "ast.all")
head(ast.all)

olat.all = ReadMotifFile(paste0("Motif_enrichment_results/Nfur_Master_Up/", overlap, "/" , peaks, "/Olat_Neo_Combine/knownResults_processed.csv"), "olat.all")
head(olat.all)

drer.all = ReadMotifFile(paste0("Motif_enrichment_results/Nfur_Master_Up/", overlap, "/" , peaks, "/Drer_Neo_Combine/knownResults_processed.csv"), "drer.all")
head(drer.all)

# Combine them into one file by finalGene column
merged_all_Up <- Reduce(function(x, y) merge(x, y, by="finalGene", all = TRUE), list(nfur.all, aaul.all, ast.all, olat.all, drer.all))
head(merged_all_Up)

merged_all_Up  = as.data.frame(merged_all_Up )
rownames(merged_all_Up) = merged_all_Up$finalGene

# Select the ones that are significant at FDR 0.1 in N. furzeri
toPlot = merged_all_Up[(merged_all_Up$nfur.all.fdr_1 <= 0.1),]
head(toPlot)
colnames(head(toPlot))

# Now get 3 sets of genes nfur only, killifish shared and conserved in Olat Drer motifs ----------------------------------------
nfur = toPlot[(toPlot$aaul.all.fdr_1 > 0.1 &
                 toPlot$ast.all.fdr_1 > 0.1 &
                 toPlot$olat.all.fdr_1 > 0.1 &
                 toPlot$drer.all.fdr_1 > 0.1
),]

killi = toPlot[((toPlot$aaul.all.fdr_1 <= 0.1 | toPlot$ast.all.fdr_1 <= 0.1) &
                  toPlot$olat.all.fdr_1 > 0.1 & toPlot$drer.all.fdr_1 > 0.1
),]

conserved = toPlot[(toPlot$olat.all.fdr_1 <= 0.1 | toPlot$drer.all.fdr_1 <= 0.1
),]

# Merge and save in anew file
toPlotNew = rbind(nfur, killi, conserved)

# Below condition must be true
sum (length(nfur[,1]), length(killi[,1]), length(conserved[,1])) == length(toPlot[,1])

write.csv(toPlotNew, file = paste0("Output_files_merged/", preFix, ".csv"))

# Mark only the top motifs to be labelled on the plot. These are selected significant labels that we want to highlight
genesLabels = toPlotNew$finalGene
index = match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARG(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels)
match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARA(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels) # Nothing should be NA
custom_label_genes = rowAnnotation(foo = anno_mark(at = index,
                                                   side = "left", # To put labels on left or riht side
                                   labels = genesLabels[index]))

# Print file with labels
pdf(paste0("Output_plots/", preFix, "_selectedLabs.pdf"), width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11, 19, 27, 35)]), 
        col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "grey80")), 
        cluster_rows = F, 
        cluster_columns = F,
        show_row_names = F,
        left_annotation = custom_label_genes # To put labels on left or riht side
)
dev.off()

# Plot another heatmap with all labels 
pdf(paste0("Output_plots/", preFix, "_AllLabs.pdf"), width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11, 19, 27, 35)]),
        col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "grey80")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 7),
        row_names_max_width = unit(30, "cm"))
dev.off()
