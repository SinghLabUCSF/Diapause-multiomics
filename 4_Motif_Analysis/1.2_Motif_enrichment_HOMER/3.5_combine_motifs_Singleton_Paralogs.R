# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ComplexHeatmap) # For Heatmap function
library(circlize) # for generating colors using colorremp2
library(gplots) # for heatmap.2 which has better row names
library(dplyr)

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

rm (list = setdiff(ls(), "ReadMotifFile"))

# Read motif enrichment files ------------------------------------------------------------------------------------------------------------------
# Global Nfur 
nfur.neoF = ReadMotifFile("Motif_enrichment_results/Nfur_Master_Up/Nfur_Neo_Combine/knownResults_processed.csv", "Paralogs")
head(nfur.neoF)

# Singletons Nfur
nfur.para = ReadMotifFile("Motif_enrichment_results/Combined_DE_Singletons/knownResults_processed.csv", "Singleton")
head(nfur.para)

merged <- Reduce(function(x, y) merge(x, y, by="finalGene", all = TRUE), list(nfur.neoF, nfur.para))
head(merged)

merged  = as.data.frame(merged)
rownames(merged) = make.unique(merged$finalGene)
head(merged)

filter1 = merged[(merged$Paralogs.fdr_1 <= 0.1 &  merged$Singleton.fdr_1 > 0.1 
                  ),]

filter2 = merged[((merged$Paralogs.fdr_1 > 0.1 & merged$Singleton.fdr_1 <= 0.1)
                  ),]

filter3 = merged[((merged$Paralogs.fdr_1 <= 0.1 & merged$Singleton.fdr_1 <= 0.1)
),]

toPlotNew = rbind(filter1, filter2, filter3)
write.csv(toPlotNew, file = "Output_files_merged/Nfur-Single-Para_Combined.csv")

# Mark only the top motifs to be labelled on the plot. These are selected significant labels that we want to highlight
genesLabels = toPlotNew$finalGene
index = match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARG(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels)
match(c("REST","FOXA1", "FOXO3", "NR2F2(NR)", "TEAD2", "HNF4A(NR)", "PPARA(NR)", "PPARA(NR)", "JUN", "CTCF", "NR2C2(NR)", "BACH2", "HRE", "KLF5"), genesLabels) # Nothing should be NA
custom_label_genes = rowAnnotation(foo = anno_mark(at = index,
                                                   side = "right", # To put labels on left or riht side
                                                   labels = genesLabels[index]))

toPlotNew = toPlotNew[complete.cases(toPlotNew),]

# Print file with labels
pdf("Output_plots/Nfur_Singleton-Para_selectedLabs.pdf", width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11)]), 
  col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "grey80")), 
  cluster_rows = F, 
  cluster_columns = F,
  show_row_names = F,
  right_annotation = custom_label_genes # To put labels on left or riht side
)
dev.off()

# Plot another heatmap with all labels 

pdf("Output_plots/Nfur_Singleton-Para_allLabs.pdf", width = 4, height = 15)
Heatmap(as.matrix(toPlotNew[, c(3, 11)]),
        col = colorRamp2(c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.45), c("orangered", "tomato", "coral", "orange", "bisque", "gray90", "gray80")), 
        #col = colorRamp2(c(0, 0.01, 0.05, 0.1, 0.25), c("orangered", "tomato", "orange", "grey80", "gray90")), 
        cluster_rows = F, 
        cluster_columns = T,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(30, "cm"))
dev.off()


