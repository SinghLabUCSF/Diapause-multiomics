# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(dplyr) # for data manipulation
library(ggfortify) # for autoplot function to plot pca
library(matrixTests) # to do the test row wise
library(qvalue)

# Protein concentrations
prot = read.table("data/embryo_protein_concentration.txt", header = T, sep = "\t")
head(prot)

# Read the table with all the conditions and tissues, modified from Yiwen's file
fullTable = read.table(file = "output/embryo_lipids_classes.txt", header = T, row.names = 1, sep = "\t")

colnames(fullTable)
colnames(fullTable) == prot$Id
colnames(fullTable) = prot$sample.name
colnames(fullTable)

# Omit one development sample that behaves weirdly in PCA
drop <- c("nfur.NonD.1")
fullTable = fullTable[ , !(names(fullTable) %in% drop)]
colnames(fullTable)

# FC is with respect to the first sample
Ttest_welch("nfur.D1m", "NonD", fullTable)
Ttest_welch("nfur.D1m", "nfur.PreD.Y", fullTable)
Ttest_welch("nfur.D1m", "nfur.Dexit", fullTable)
Ttest_welch("nfur.D1m", "nfur.D6d", fullTable)
Ttest_welch("nfur.D6d", "NonD", fullTable)
Ttest_welch("nfur.D6d", "nfur.PreD.Y", fullTable)
Ttest_welch("nfur.D6d", "nfur.Dexit", fullTable)
Ttest_welch("NonD", "nfur.PreD.Y", fullTable)
Ttest_welch("NonD", "nfur.Dexit", fullTable)
Ttest_welch("nfur.PreD.Y", "Ast.PreD.Y", fullTable)
Ttest_welch("nfur.PreD.Y", "nfur.PreD.O", fullTable)
Ttest_welch("nfur.PreD.Y", "nfur.Dexit", fullTable)

# Input to this function are cond1 and cond2 and fullTable
Ttest_welch <- function(cond1, cond2, Lipids) {
  
  c1 = dplyr::select(Lipids, matches(cond1))
  c2 = dplyr::select(Lipids, matches(cond2))
  print (head(c1))
  print (head(c2))
  
  c1_c2 = row_t_welch(c1, c2)
  c1_c2$padj = p.adjust(c1_c2$pvalue, method = "BH")
  
  #qvalues = qvalue(c1_c2$pvalue)
  #c1_c2$qval = qvalues$qvalues
  
  c1_c2$FC = log2(c1_c2$mean.x / c1_c2$mean.y)
  
  final = data.frame(class = rownames(Lipids), c1, c2, pvalue = c1_c2$pvalue, padj = c1_c2$padj, FC = c1_c2$FC)
  write.table(final, file = paste0("DE_lipid_classes/", cond1, "_v_", cond2, ".txt"), row.names = F, quote = F, sep = "\t")
  
  rm(c1, c2, cond1, cond2, c1_c2, final)
}

