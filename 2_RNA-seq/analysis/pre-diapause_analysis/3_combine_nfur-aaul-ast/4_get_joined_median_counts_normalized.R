# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read normalized tpm values
tpm = read.table("CountsNormDESeq2_aaul-ast-nfur.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)
colnames(tpm)
# Combined developmental stages--see readme
colnames(dplyr::select(tpm, matches("nfur_kv_Y|Y6_KV_R2")))
tpm = mutate(tpm, median.kvY.nfur = rowMedians(as.matrix(dplyr::select(tpm, matches("nfur_kv_Y|Y6_KV_R2"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("Aaul_kv")))
tpm = mutate(tpm, median.kv.aaul = rowMedians(as.matrix(dplyr::select(tpm, matches("Aaul_kv"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("Ast_kv_Y")))
tpm = mutate(tpm, median.kvY.ast = rowMedians(as.matrix(dplyr::select(tpm, matches("Ast_kv_Y"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

write.table(tpm, file = "Normalized-Median-Counts_nfur-aaul-ast.csv", row.names = F, sep = ",")

