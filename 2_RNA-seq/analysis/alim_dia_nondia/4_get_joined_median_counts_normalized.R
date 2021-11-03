# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read normalized tpm values
#tpm = read.table("CountsNormDESeq2_alim_dia-nondia.csv",header = T, sep = ",")
tpm = read.table("Counts_alim_dia-nondia.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)
colnames(tpm)
# Combined developmental stages--see readme
colnames(dplyr::select(tpm, matches("alim.dia")))
tpm = mutate(tpm, median.dia.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.dia"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("4dpdD")))
tpm = mutate(tpm, median.4dpdD.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("4dpdD"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

#write.table(tpm, file = "Normalized-Median-Counts_alim_dia-nondia.csv", row.names = F, sep = ",")
write.table(tpm, file = "Raw-Median-Counts_alim_dia-nondia.csv", row.names = F, sep = ",")
