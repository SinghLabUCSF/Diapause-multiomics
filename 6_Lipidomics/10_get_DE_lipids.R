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

# input file Phospholipids or AllLipids or Protein
input = "Phospholipids"

# Read the table with all the conditions
fullTable = read.table(file = paste0("output/embryo_lipids_normalized_", input, ".txt"), header = T, row.names = 1, sep = "\t")

colnames(fullTable)
colnames(fullTable)[5:36] == prot$Id
colnames(fullTable)[5:36] = prot$sample.name
colnames(fullTable)

# Omit one development sample that behaves weirdly in PCA
drop <- c("nfur.NonD.1")
fullTable = fullTable[ , !(names(fullTable) %in% drop)]
colnames(fullTable)

Ttest_welch("nfur.D1m", "NonD", fullTable)
Ttest_welch("nfur.D1m", "nfur.PreD.Y", fullTable)
Ttest_welch("nfur.D1m", "nfur.Dexit", fullTable)
Ttest_welch("nfur.D1m", "nfur.D6d", fullTable)
Ttest_welch("nfur.D6d", "NonD", fullTable)
Ttest_welch("nfur.D6D", "nfur.PreD.Y", fullTable)
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
  
  final = data.frame(m.z = Lipids$m.z, class = Lipids$Class, FattyAcid = Lipids$FattyAcid, c1, c2, pvalue = c1_c2$pvalue, padj = c1_c2$padj, FC = c1_c2$FC)
  write.table(final, file = paste0("DE_individual_lipids/", cond1, "_v_", cond2, ".txt"), row.names = F, quote = F, sep = "\t")
  
  rm(c1, c2, cond1, cond2, c1_c2, final)
}

