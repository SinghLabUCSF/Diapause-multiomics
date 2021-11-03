# Join differential expession lists
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read the DE files
difDia = read.table(file="DE_Dia-vs-4dpD.txt", head = T)

# Process and add gene columns
difDia = cbind(genes = rownames(difDia), difDia)
difDia = difDia[, c(1,3,7)]
colnames(difDia) = c("genes", "Dia.log2FC", "Dia.padj")
head(difDia)


# Read normalized tpm values
tpm = read.table("CountsNormDESeq2_alim_dia-nondia.csv",header = T, sep = ",")
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)

colnames(dplyr::select(tpm, matches("4dpd")))
tpm = mutate(tpm, median.Development.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("4dpd"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("alim.dia")))
tpm = mutate(tpm, median.Diapause.alim = rowMedians(as.matrix(dplyr::select(tpm, matches("alim.dia"))), na.rm = TRUE))

# Only take median columns -- I am keeping all the columns for now
#tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

finalMerged <- Reduce(function(x, y) merge(x, y, by="genes"), list(difDia, tpm))
colnames(finalMerged)
head(finalMerged)

write.table(finalMerged, file = "Counts_DE_Combined_Final_alim_dia-nondia.csv", row.names = F, sep = ",")

# Now make the lists of up and down genes
######################## UP GENES #########################################
# List 2: Upregulated wrt both preD and nonD with FDR 0.01 and FC > 1
dpUp1 = finalMerged[(finalMerged$Dia.log2FC >= 0 & finalMerged$Dia.padj < 0.05),]
length(dpUp1[,1])  # 7362
write.table(dpUp1, file = "dpUp.txt", row.names = F, sep = "\t")

################### DOWN ###############################
# List 2: Significantly down wrt both preD and nonD with FDR 0.05
dpDown1 = finalMerged[(finalMerged$Dia.log2FC < 0 & finalMerged$Dia.padj < 0.05),]
length(dpDown1[,1]) # 5320
write.table(dpDown1, file = "dpDown.txt", row.names = F, sep = "\t")


# Plot boxplots for most relaxed data
boxplot(dpUp1$Dia.log2FC, dpDown1$Dia.log2FC, outline = F)
boxplot(dpUp1$median.Diapause.alim, dpUp1$median.Development.alim, dpDown1$median.Diapause.alim, dpDown1$median.Development.alim, outline = F)

boxplot(dpUp3$median.Diapause, dpUp3$median.Development, outline = F)
boxplot(dpDown3$median.Diapause, dpDown3$median.Development, outline = F)
