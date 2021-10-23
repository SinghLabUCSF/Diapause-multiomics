library("matrixStats")

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

lipids = read.table("embryo_lipids_normalized_IS.txt", header = T, sep = "\t")
head(lipids)

# Check which Phospholipids are there (PA, PE, PC, PI, PS, PG)
unique(lipids[,43])

# Calculate normalization factor based only on Phospholipids
medians = colMedians(as.matrix(
  lipids[lipids$Class %in% c("PA", "PE", "PC", "PI", "PS", "PG"), 6:37]
  )) 

col_medians = as.data.frame(medians)
rownames(col_medians) = colnames(lipids[,6:37])
col_medians
write.table(col_medians, "output/embryo_medinas_after_normalization_Phospholipids.txt",quote = F, row.names = T, col.names = T, sep = "\t")


global_median = median(medians)
normfactors = medians/global_median
normfactors
length(normfactors)

normfactorsToWrite = as.data.frame(normfactors)
row.names(normfactorsToWrite) = colnames(lipids[,6:37])
write.table(normfactorsToWrite, "output/embryo_normalization_factors_Phospholipids.txt",quote = F, row.names = T, col.names = T, sep = "\t")

# To test if mapply is working fine here
#mapply(`*`, lipids[1:5,6:40], 0:34)

norm_vals = mapply(`/`, lipids[,6:37], normfactors)


# ### Xiaoai Code-checking ===================================================================================
# #normalize by dividing raw data with normalization factor, so that the new column median in the normalized matrix should be identical
# 
# #subset matrix to get phospholipids only
# mtx_phos <- as.matrix(
#   +     lipids[lipids$Class %in% c("PA", "PE", "PC", "PI", "PS", "PG"), 6:37])
# #Raw matrix divided by normalization factor to perform normalization
# norm_phos.val <- mapply('/', as.data.frame(mtx_phos), normfactors)# divide instead of multiply
# #After normalization, column median of the new matrix should be identical across all samples
# phos.median <- colMedians(norm_phos.val)
# phos.median
# ### ========================================================================================================


# Bind and make the final file with normalized values
normalized_lipids = cbind(lipids[,1:5], norm_vals, lipids[,38:48])
head(normalized_lipids)
write.table(normalized_lipids, "output/embryo_lipids_normalized_Phospholipids.txt",quote = F, row.names = F, col.names = T, sep = "\t")

# Change the headers to more meaningful values from protein file. This is for people to have a look
prot = read.table("embryo_protein_concentration.txt", header = T, sep = "\t")
selected_cols = normalized_lipids[, 6:37] # Take lipid columns
colnames(selected_cols) == prot$Id # Check if the order is same
colnames(selected_cols) = prot$sample.name # Add prper names
head(selected_cols) # Check
# Add the 3 lipid ion columns
selected_cols$LipidIon = normalized_lipids$LipidIon
selected_cols$FattyAcid = normalized_lipids$FattyAcid
selected_cols$Class = normalized_lipids$Class
head(selected_cols)
colnames(selected_cols)

# Change the order and select columns
selected_cols = selected_cols[, c("LipidIon", "FattyAcid", "Class", 
                                  "nfur.PreD.Y.1", "nfur.PreD.Y.2", "nfur.PreD.Y.3", "nfur.PreD.Y.4",
                                  "nfur.NonD.2", "nfur.NonD.3", "nfur.NonD.4",
                                  "nfur.D6d.1", "nfur.D6d.2", "nfur.D6d.3", "nfur.D6d.4", 
                                  "nfur.D1m.1", "nfur.D1m.2", "nfur.D1m.3", "nfur.D1m.4", 
                                  "Ast.PreD.Y.1", "Ast.PreD.Y.2", "Ast.PreD.Y.3", "Ast.PreD.Y.4",
                                  "nfur.Dexit.1", "nfur.Dexit.2", "nfur.Dexit.3", "nfur.Dexit.4"
                                  )] 
write.table(selected_cols, "output/embryo_lipids_normalized_Phospholipids_WithProperNames.txt",quote = T, row.names = F, col.names = T, sep = "\t")

