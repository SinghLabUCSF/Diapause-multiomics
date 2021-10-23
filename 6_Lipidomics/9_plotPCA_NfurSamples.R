# To generate all the PCA combinations for lipidomics data
# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr) # for data manipulation
library(ggfortify) # for autoplot function to plot pca

# Protein concentrations
prot = read.table("data/embryo_protein_concentration.txt", header = T, sep = "\t")
head(prot)

# input file Phospholipids or AllLipids or Protein
input = "Phospholipids"

# Read the table with all the conditions
fullTable = read.table(file = paste0("output/embryo_lipids_normalized_", input, ".txt"), header = T, row.names = 1, sep = "\t")
fullTable = fullTable[apply(fullTable!=0, 1, all),] # 0 values can create issues for Ast samples so I remove any rows containg 0

fullTable = fullTable[, 5:36] # The column number is one less because first one is rowname
colnames(fullTable)
colnames(fullTable) == prot$Id
colnames(fullTable) = prot$sample.name

# Now filter for different criteria and make PCA plots
sampleData = dplyr::select(fullTable, matches("nfur.D1m|nfur.D6d|NonD.2|NonD.3|NonD.4|Nfur.PreD.Y"))

colnames(sampleData) # check columns
trMatrix = t(na.omit(sampleData)) # transpose

rownames(trMatrix)
colnames(trMatrix)

# Add some metadata to change colors etc.
finalTable = cbind.data.frame(sample.name  = rownames(trMatrix),
                              sampleId  = prot$Id[prot$sample.name %in% rownames(trMatrix)],
                              species  = prot$species[prot$sample.name %in% rownames(trMatrix)],
                              stage  = prot$stage[prot$sample.name %in% rownames(trMatrix)],
                              trMatrix)
#forPCA = finalTable
finalTable[, 1:5] # To check how many columns need to be excluded
forPCA = finalTable[,5:length(names(finalTable))]
pca <- prcomp(forPCA, center = TRUE, scale. = T)
summary(pca)

pdf (file = "PCA_Nfur-Samples.pdf", width = 6, height = 4)
autoplot(pca,
         x = 1, y = 2, # For other PCs: e.g. PC2, PC3
         data = finalTable,
         colour = 'stage',
         label = TRUE,
         shape = 'species',
         #frame = TRUE,
         #frame.type = 'norm',
         size = 6
         #alpha=I(0.5)
) + theme_bw() 
dev.off()
