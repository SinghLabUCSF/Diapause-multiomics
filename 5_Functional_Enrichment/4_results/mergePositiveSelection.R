# To merge the lipis related functions for RNA and ATAC
# Input are functional enrichment outputs
# IMPORTANT: From the output file remove ' and relcace by nothing, else the numbers will be wrong. I do it manually.
 
# Set wd to the current directory
library(ComplexHeatmap)
library(circlize) # for generating colors
library(gplots)
library(pheatmap)
library (dplyr)
library("data.table") # It has fread function that can take care of
require(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read ATAC files
goBP_atac = fread("GOBP_PositiveSelection_Peaks.txt", sep = "\t", header = T)
goMF_atac = fread("GOMF_PositiveSelection_Peaks.txt", sep = "\t", header = T)
goCC_atac = fread("GOCC_PositiveSelection_Peaks.txt", sep = "\t", header = T)
print ("*** Check the length of file is okay before proceeding. ***")
# If length differ than the terms in the list, remove ' or other special characters from the Go name and rerun

# Select only relevant columns
goBP_atac = goBP_atac[, c(1,2,7)]
goMF_atac = goMF_atac[, c(1,2,7)]
goCC_atac = goCC_atac[, c(1,2,7)]

# Replace first colum name
colnames(goBP_atac) = c("GOID", "BP.Pvalue.ATAC", "BP.Term.ATAC")
goBP_atac$Ontology = "BP" # Add an ontology column

colnames(goMF_atac) = c("GOID", "MF.Pvalue.ATAC", "MF.Term.ATAC")
goMF_atac$Ontology = "MF" # Add an ontology column

colnames(goCC_atac) = c("GOID", "CC.Pvalue.ATAC", "CC.Term.ATAC")
goCC_atac$Ontology = "CC" # Add an ontology column

head(goBP_atac)
head(goMF_atac)
head(goCC_atac)

# Merge ATAC terms ----------------------------------------
Atac <- Reduce(function(x, y) merge(x, y, by="GOID", all = TRUE), list(goBP_atac, goMF_atac, goCC_atac))
head(Atac)

Ontology.ATAC = coalesce(Atac$Ontology.x, Atac$Ontology.y, Atac$Ontology)
head(Ontology.ATAC)

Pvalue.ATAC = coalesce(Atac$BP.Pvalue.ATAC, Atac$MF.Pvalue.ATAC, Atac$CC.Pvalue.ATAC)
head(Pvalue.ATAC)

Term.ATAC = coalesce(Atac$BP.Term.ATAC, Atac$MF.Term.ATAC, Atac$CC.Term.ATAC)
head(Term.ATAC)

head(Atac)

Atac$Ontology.ATAC = Ontology.ATAC
Atac$Pvalue.ATAC = Pvalue.ATAC
Atac$Term.ATAC = Term.ATAC
head(Atac)

Atac_merged = Atac[, c("GOID", "Ontology.ATAC" , "Pvalue.ATAC", "Term.ATAC")]
head(Atac_merged)

lipidLogical = Atac_merged$Term %like% "lipid|fatty|glycer| fat |lipase"
Atac_merged$lipid.related = lipidLogical
head(Atac_merged)

write.csv(Atac_merged, file = "GOALL_ATAC_PositiveSelection.csv")
rownames(Atac_merged) = make.unique(Atac_merged$Term.ATAC)

Lipid_Atac_merged = Atac_merged[Atac_merged$Term.ATAC %like% "lipid|fatty|glycer| fat |lipase",]
head(Lipid_Atac_merged)

toPlot = Lipid_Atac_merged[which(Lipid_Atac_merged$Pvalue.ATAC <= 0.05), ]
head(toPlot)
write.csv(toPlot, file = "GOALL_ATAC_PositiveSelection_LipidRelated.csv")

toPlotNew = Atac_merged[which(Atac_merged$Pvalue.ATAC <= 0.05), ]
head(toPlotNew)
write.csv(toPlotNew, file = "GOALL_ATAC_PositiveSelection_Significant.csv")

