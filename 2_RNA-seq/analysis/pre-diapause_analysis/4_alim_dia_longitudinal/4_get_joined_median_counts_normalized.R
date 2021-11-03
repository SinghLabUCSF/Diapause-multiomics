# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read normalized tpm values
tpm = read.table("CountsNormDESeq2_alim_dia-longitudianl.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)
colnames(tpm)

# ---------------------- for 30C embryos ----------------------------------
colnames(dplyr::select(tpm, matches("alim.30C_Neural_keel")))
tpm = mutate(tpm, median.30C_Neural_keel.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_Neural_keel"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_dispersed")))
tpm = mutate(tpm, median.30C_dispersed.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_dispersed"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_6_somites")))
tpm = mutate(tpm, median.30C_6_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_6_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_10_somites")))
tpm = mutate(tpm, median.30C_10_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_10_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_16_somites")))
tpm = mutate(tpm, median.30C_16_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_16_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_20_somites")))
tpm = mutate(tpm, median.30C_20_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_20_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_24_somites")))
tpm = mutate(tpm, median.30C_24_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_24_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.30C_6_somites|alim.30C_10_somites|alim.30C_16_somites|alim.30C_20_somites|alim.30C_24_somites")))
tpm = mutate(tpm, median.30C_AllSomites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.30C_6_somites|alim.30C_10_somites|alim.30C_16_somites|alim.30C_20_somites|alim.30C_24_somites"))), na.rm = TRUE))

# ---------------------- for 20C embryos ----------------------------------
colnames(dplyr::select(tpm, matches("alim.20C_Neural_keel")))
tpm = mutate(tpm, median.20C_Neural_keel.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_Neural_keel"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_dispersed")))
tpm = mutate(tpm, median.20C_dispersed.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_dispersed"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_6_somites")))
tpm = mutate(tpm, median.20C_6_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_6_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_10_somites")))
tpm = mutate(tpm, median.20C_10_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_10_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_16_somites")))
tpm = mutate(tpm, median.20C_16_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_16_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_20_somites")))
tpm = mutate(tpm, median.20C_20_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_20_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_24_somites")))
tpm = mutate(tpm, median.20C_24_Somites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_24_somites"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.20C_6_somites|alim.20C_10_somites|alim.20C_16_somites|alim.20C_20_somites|alim.20C_24_somites")))
tpm = mutate(tpm, median.20C_AllSomites.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.20C_6_somites|alim.20C_10_somites|alim.20C_16_somites|alim.20C_20_somites|alim.20C_24_somites"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

write.table(tpm, file = "Normalized-Median-Counts_alim-longitudinal.csv", row.names = F, sep = ",")

