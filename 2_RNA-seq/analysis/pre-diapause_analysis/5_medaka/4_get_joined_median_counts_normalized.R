# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read normalized tpm values
tpm = read.table("CountsNormDESeq2_medaka_development_all.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)

colnames(dplyr::select(tpm, matches("St11")))
tpm = mutate(tpm, median.St11.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St11"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("St13")))
tpm = mutate(tpm, median.St13.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St13"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("St16")))
tpm = mutate(tpm, median.St16.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St16"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("St19")))
tpm = mutate(tpm, median.St19.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St19"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("St25")))
tpm = mutate(tpm, median.St25.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St25"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("St32")))
tpm = mutate(tpm, median.St32.olat = rowMedians(as.matrix(dplyr::select(tpm, matches("St32"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

write.table(tpm, file = "Normalized-Median-Counts_medaka_all.csv", row.names = F, sep = ",")

