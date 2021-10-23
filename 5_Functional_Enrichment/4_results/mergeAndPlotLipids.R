# To merge the lipis related functions for RNA and ATAC
# Input are functional enrichment outputs
 
# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ComplexHeatmap)
library(circlize) # for generating colors
library(gplots)
library(pheatmap)
library (dplyr)
library("data.table") 
require(gridExtra)

# Read RNA files
goBP_neof1 = fread("GOBP_NeoF1_Genes.txt", sep = "\t", header = T)
goMF_neof1 = fread("GOMF_NeoF1_Genes.txt", sep = "\t", header = T)
goCC_neof1 = fread("GOCC_NeoF1_Genes.txt", sep = "\t", header = T)

print ("*** Check the length of file is okay before proceeding. ***")

# Read ATAC files
goBP_atac = fread("GOBP_DiapauseUpAll_Peaks.txt", sep = "\t", header = T)
goMF_atac = fread("GOMF_DiapauseUpAll_Peaks.txt", sep = "\t", header = T)
goCC_atac = fread("GOCC_DiapauseUpAll_Peaks.txt", sep = "\t", header = T)

print ("*** Check the length of file is okay before proceeding. ***")
# If length differ than the terms in the list, remove ' or other special characters from the GO name and rerun

# Select only relevant columns
goBP_neof1 = goBP_neof1[, c(1,2,7)]
goMF_neof1 = goMF_neof1[, c(1,2,7)]
goCC_neof1 = goCC_neof1[, c(1,2,7)]

goBP_atac = goBP_atac[, c(1,2,7)]
goMF_atac = goMF_atac[, c(1,2,7)]
goCC_atac = goCC_atac[, c(1,2,7)]

# Replace first column name
colnames(goBP_neof1) = c("GOID", "BP.Pvalue.RNA", "BP.Term.RNA")
goBP_neof1$Ontology = "BP" # Add an ontology column

colnames(goMF_neof1) = c("GOID", "MF.Pvalue.RNA", "MF.Term.RNA")
goMF_neof1$Ontology = "MF" # Add an ontology column

colnames(goCC_neof1) = c("GOID", "CC.Pvalue.RNA", "CC.Term.RNA")
goCC_neof1$Ontology = "CC" # Add an ontology column

colnames(goBP_atac) = c("GOID", "BP.Pvalue.ATAC", "BP.Term.ATAC")
goBP_atac$Ontology = "BP" # Add an ontology column

colnames(goMF_atac) = c("GOID", "MF.Pvalue.ATAC", "MF.Term.ATAC")
goMF_atac$Ontology = "MF" # Add an ontology column

colnames(goCC_atac) = c("GOID", "CC.Pvalue.ATAC", "CC.Term.ATAC")
goCC_atac$Ontology = "CC" # Add an ontology column

head(goBP_neof1)
head(goMF_neof1)
head(goCC_neof1)

head(goBP_atac)
head(goMF_atac)
head(goCC_atac)

# Merge NeoF1 terms ----------------------------------------
neoF1 <- Reduce(function(x, y) merge(x, y, by="GOID", all = TRUE), list(goBP_neof1, goMF_neof1, goCC_neof1))
head(neoF1)

Ontology.RNA = coalesce(neoF1$Ontology.x, neoF1$Ontology.y, neoF1$Ontology)
head(Ontology.RNA)

Pvalue.RNA = coalesce(neoF1$BP.Pvalue.RNA, neoF1$MF.Pvalue.RNA, neoF1$CC.Pvalue.RNA)
head(Pvalue.RNA)

Term.RNA = coalesce(neoF1$BP.Term.RNA, neoF1$MF.Term.RNA, neoF1$CC.Term.RNA)
head(Term.RNA)

head(neoF1)

neoF1$Ontology.RNA = Ontology.RNA
neoF1$Pvalue.RNA = Pvalue.RNA
neoF1$Term.RNA = Term.RNA
head(neoF1)

neoF1_merged = neoF1[, c("GOID", "Ontology.RNA" , "Pvalue.RNA", "Term.RNA")]
head(neoF1_merged)


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

# Merge RNA and ATAC and reduce the columns again
RNA_ATAC_merge <- Reduce(function(x, y) merge(x, y, by="GOID", all = TRUE), list(neoF1_merged, Atac_merged))
head(RNA_ATAC_merge)

Ontology = coalesce(RNA_ATAC_merge$Ontology.RNA, RNA_ATAC_merge$Ontology.ATAC)
head(Ontology)

Term = coalesce(RNA_ATAC_merge$Term.RNA, RNA_ATAC_merge$Term.ATAC)
head(Term)

RNA_ATAC_merge$Ontology = Ontology
RNA_ATAC_merge$Term = Term

head(RNA_ATAC_merge)

RNA_ATAC_merge = RNA_ATAC_merge[, c("GOID", "Term", "Pvalue.RNA", "Pvalue.ATAC", "Ontology")]
head(RNA_ATAC_merge)

lipidLogical = RNA_ATAC_merge$Term %like% "lipid|fatty|glycer| fat |lipase"
RNA_ATAC_merge$lipid.related = lipidLogical
head(RNA_ATAC_merge)

write.csv(RNA_ATAC_merge, file = "GOALL_RNA_ATAC_merged.csv")
rownames(RNA_ATAC_merge) = make.unique(RNA_ATAC_merge$Term)
colnames(RNA_ATAC_merge)

Lipid_RNA_ATAC_merge = RNA_ATAC_merge[RNA_ATAC_merge$Term %like% "lipid|fatty|glycer| fat |lipase",] 
head(Lipid_RNA_ATAC_merge)

toPlot = Lipid_RNA_ATAC_merge[which(Lipid_RNA_ATAC_merge$Pvalue.RNA <= 0.05 & Lipid_RNA_ATAC_merge$Pvalue.ATAC <= 0.05), ]
head(toPlot)
write.csv(toPlot, file = "GOALL_RNA_ATAC_merged_LipidRelated.csv")

toPlotNew = RNA_ATAC_merge[which(RNA_ATAC_merge$Pvalue.RNA <= 0.05 & RNA_ATAC_merge$Pvalue.ATAC <= 0.05), ]
head(toPlotNew)
write.csv(toPlotNew, file = "GOALL_RNA_ATAC_merged_AllSignifiant.csv")

pdf(file = "RNA_ATAC_lipid_bar.pdf", height = 2, width = 10)
plot1 = ggplot(data = toPlot, aes(x = reorder(Term, -Pvalue.ATAC), y = -log10(Pvalue.ATAC))) +
        geom_bar(stat = 'identity') +
        coord_flip() +
        theme_bw()  

plot2 = ggplot(data = toPlot, aes(x = reorder(Term, -Pvalue.ATAC), y = -log10(Pvalue.RNA))) +
        geom_bar(stat = 'identity') +
        coord_flip() +
        theme_bw() 
#        theme(axis.text.y=element_blank()) +
#        theme(axis.title.y=element_blank())

grid.arrange(plot1, plot2, ncol=2)
dev.off()

# plot3 = ggplot(data = toPlot, aes(y = reorder(Term, -Pvalue.ATAC), x = -log10(Pvalue.ATAC)), size = -log10(Pvalue.ATAC)) +
#         geom_point(stat = 'identity') +
#         scale_size(range = c(1, 10), name="Pvaue") +
#         theme_bw() 
# 
# plot4 = ggplot(data = toPlot, aes(y = reorder(Term, -Pvalue.ATAC), x = -log10(Pvalue.RNA))) +
#         geom_point(stat = 'identity') +
#         theme_bw() 
# 
# grid.arrange(plot3, plot4, ncol=2)

# Heatmap(as.matrix(toPlot[,c(3,4)]),
#         cluster_rows = F, 
#         cluster_columns = F,
#         show_row_names = T,
#         row_names_max_width = unit(10, "cm"))

