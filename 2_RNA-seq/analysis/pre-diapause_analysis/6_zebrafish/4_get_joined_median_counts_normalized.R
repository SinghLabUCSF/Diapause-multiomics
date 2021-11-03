# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read normalized tpm values
tpm = read.table("CountsNormDESeq2_zebrafish_development_all.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)
colnames(tpm)
# Combined developmental stages--see readme
colnames(dplyr::select(tpm, matches("2.4cell|1Kcell|X2_hpf")))
tpm = mutate(tpm, median.0to3hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("2.4cell|1Kcell|X2_hpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("dome|shield|6hpf|X6_hpf|X8_hpf")))
tpm = mutate(tpm, median.4to8hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("dome|shield|6hpf|X6_hpf|X8_hpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("bud|X16_hpf|X12_hpf")))
tpm = mutate(tpm, median.10to16hpf.mid.somites.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("bud|X16_hpf|X12_hpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("X20_hpf|X26_hpf")))
tpm = mutate(tpm, median.20to26hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("X20_hpf|X26_hpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("X28hpf|X30_hpf")))
tpm = mutate(tpm, median.28to30hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("X28hpf|X30_hpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("X36_hpf|X48_hpf|X72_hpf|X2dpf")))
tpm = mutate(tpm, median.36to72hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("X36_hpf|X48_hpf|X72_hpf|X2dpf"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("X5dpf")))
tpm = mutate(tpm, median.120hpf.drer = rowMedians(as.matrix(dplyr::select(tpm, matches("X5dpf"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

write.table(tpm, file = "Normalized-Median-Counts_zebrafish_all.csv", row.names = F, sep = ",")

