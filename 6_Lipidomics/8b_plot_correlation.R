library(ggplot2)
library(ggrepel)

# Set wd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Protein concentrations
prot = read.table("data/embryo_protein_concentration.txt", header = T, sep = "\t")
head(prot)

# All lipids
# allLipids = read.table("embryo_medinas_after_normalization_AllLipids.txt", header = T, sep = "\t")
# allLipids$Id = rownames(allLipids)
# colnames(allLipids) = c("median.allLipids", "Id")
# head(allLipids)

# Phospholipids
phosphoLipids = read.table("output/embryo_medinas_after_normalization_Phospholipids.txt", header = T, sep = "\t")
phosphoLipids$Id = rownames(phosphoLipids)
colnames(phosphoLipids) = c("median.phosphoLipids", "Id")
head(phosphoLipids)


# Combine the data frames
lpcor <- Reduce(function(x, y) merge(x, y, by="Id", all = TRUE), list(prot, phosphoLipids))
head(lpcor)

# Correlation between proteins and all lipids
# ggplot(lpcor, aes(x=log(median.allLipids), y=protein)) +
#   geom_point(size = 6) +
#   geom_smooth(method=lm) +
#   geom_text_repel(aes(label=lpcor$sample.name))
# 
# cor.test(log(lpcor$median.allLipids), lpcor$protein)

# Correlation between proteins and phospholipids
pdf (file = "Protein_PC_Correlation.pdf", width = 5, height = 5)
ggplot(lpcor, aes(x=log(median.phosphoLipids), y=log(protein))) + 
  geom_point(size = 6) +  
  geom_smooth(method=lm) +
  geom_text_repel(aes(label=lpcor$sample.name))
dev.off()

cor.test(log(lpcor$median.phosphoLipids), log(lpcor$protein), method = "pearson")

