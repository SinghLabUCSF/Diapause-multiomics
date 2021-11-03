# Join differential expression lists. And selectes the genes that are up and down in diapause.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse") # This includes dplyer
library("matrixStats")

# Read the pairwise DE files
difD3d = read.table(file="2_Counts_DE/DE_D3d-vs-PreD.txt", head = T)
difD6d = read.table(file="2_Counts_DE/DE_D6d-vs-PreD.txt", head = T)
difD1m = read.table(file="2_Counts_DE/DE_D1m-vs-PreD.txt", head = T)
difNonD = read.table(file="2_Counts_DE/DE_NonD-vs-PreD.txt", head = T)
difD3d2 = read.table(file="2_Counts_DE/DE_D3d-vs-nonD.txt", head = T)
difD6d2 = read.table(file="2_Counts_DE/DE_D6d-vs-nonD.txt", head = T)
difD1m2 = read.table(file="2_Counts_DE/DE_D1m-vs-nonD.txt", head = T)

# Process and add gene columns
difD3d = cbind(genes = rownames(difD3d), difD3d)
difD3d = difD3d[, c(1,3,7)]
colnames(difD3d) = c("genes", "D3d.log2FC", "D3d.padj")
head(difD3d)

difD6d = cbind(genes = rownames(difD6d), difD6d)
difD6d = difD6d[, c(1,3,7)]
colnames(difD6d) = c("genes", "D6d.log2FC", "D6d.padj")
head(difD6d)

difD1m = cbind(genes = rownames(difD1m), difD1m)
difD1m = difD1m[, c(1,3,7)]
colnames(difD1m) = c("genes", "D1m.log2FC", "D1m.padj")
head(difD1m)

difNonD = cbind(genes = rownames(difNonD), difNonD)
difNonD = difNonD[, c(1,3,7)]
colnames(difNonD) = c("genes", "nonD.log2FC", "nonD.padj")
head(difNonD)

difD3d2 = cbind(genes = rownames(difD3d2), difD3d2)
difD3d2 = difD3d2[, c(1,3,7)]
colnames(difD3d2) = c("genes", "D3d_n.log2FC", "D3d_n.padj")
head(difD3d2)

difD6d2 = cbind(genes = rownames(difD6d2), difD6d2)
difD6d2 = difD6d2[, c(1,3,7)]
colnames(difD6d2) = c("genes", "D6d_n.log2FC", "D6d_n.padj")
head(difD6d2)

difD1m2 = cbind(genes = rownames(difD1m2), difD1m2)
difD1m2 = difD1m2[, c(1,3,7)]
colnames(difD1m2) = c("genes", "D1m_n.log2FC", "D1m_n.padj")
head(difD1m2)

# Read normalized tpm values
tpm = read.table("2_Counts_DE/CountsNormDESeq2_nfur_diapause_CK.csv",header = T, sep = ",")
#tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

# Now add columns that have median values for all conditions
# Remember I am using dplyr::select because select is there is many packages and so it gives an error if thee is a conflict
# I'll print the column names that are being averaged, check that they are correct
genes = rownames(tpm)
tpm = mutate(tpm, median.overall = rowMedians(as.matrix(tpm), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("preD|nonD")))
tpm = mutate(tpm, median.Development = rowMedians(as.matrix(dplyr::select(tpm, matches("preD|nonD"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("D3d|D6d|D1m")))
tpm = mutate(tpm, median.Diapause = rowMedians(as.matrix(dplyr::select(tpm, matches("D3d|D6d|D1m"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("preD")))
tpm = mutate(tpm, median.preD = rowMedians(as.matrix(dplyr::select(tpm, matches("preD"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("D3d")))
tpm = mutate(tpm, median.D3d = rowMedians(as.matrix(dplyr::select(tpm, matches("D3d"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("D6d")))
tpm = mutate(tpm, median.D6d = rowMedians(as.matrix(dplyr::select(tpm, matches("D6d"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("D1m")))
tpm = mutate(tpm, median.D1m = rowMedians(as.matrix(dplyr::select(tpm, matches("D1m"))), na.rm = TRUE))

colnames(dplyr::select(tpm, matches("nonD")))
tpm = mutate(tpm, median.nonD = rowMedians(as.matrix(dplyr::select(tpm, matches("nonD"))), na.rm = TRUE))


# Only take median columns -- I am keeping all the columns for now
#tpm = dplyr::select(tpm, matches("median"))

colnames(tpm)

# Assign rownames again
head(tpm)
rownames(tpm) = genes
head(tpm)
tpm = cbind(genes = rownames(tpm), tpm)
head(tpm)

# Final merged file
finalMerged <- Reduce(function(x, y) merge(x, y, by="genes"), list(difD3d, difD6d, difD1m, difNonD, difD3d2, difD6d2, difD1m2, tpm))

# Median Fold Change
colnames(dplyr::select(finalMerged, matches("D3d.log2FC|D6d.log2FC|D1m.log2FC|D3d_n.log2FC|D6d_n.log2FC|D1m_n.log2FC")))
finalMerged = mutate(finalMerged, median.FC.Diapause = rowMedians(as.matrix(dplyr::select(finalMerged, matches("D3d.log2FC|D6d.log2FC|D1m.log2FC|D3d_n.log2FC|D6d_n.log2FC|D1m_n.log2FC"))), na.rm = TRUE))

# Median padj. THIS WILLBE USED TO PLOT ON THE FIGURES.
colnames(dplyr::select(finalMerged, matches("D3d.padj|D6d.padj|D1m.padj|D3d_n.padj|D6d_n.padj|D1m_n.padj")))
finalMerged = mutate(finalMerged, median.FC.padj = rowMedians(as.matrix(dplyr::select(finalMerged, matches("D3d.padj|D6d.padj|D1m.padj|D3d_n.padj|D6d_n.padj|D1m_n.padj"))), na.rm = TRUE))

colnames(finalMerged)
head(finalMerged)

write.table(finalMerged, file = "2_Counts_DE/Counts_DE_Combined_Final.csv", row.names = F, sep = ",")

# Now make two lists of up and down genes. I am trying bunch of criteria. But will use the relaxed one.
######################## UP GENES #########################################
# # List 2: Upregulated wrt both preD and nonD with FDR 0.01 and FC > 1
# dpUp1 = finalMerged[(((finalMerged$D3d.log2FC > 1 & finalMerged$D3d.padj < 0.01) &
#                         (finalMerged$D6d.log2FC > 1 & finalMerged$D6d.padj < 0.01) &
#                         (finalMerged$D1m.log2FC > 1 & finalMerged$D1m.padj < 0.01)) |
#                        ((finalMerged$D3d_n.log2FC > 1 & finalMerged$D3d_n.padj < 0.01) &
#                           (finalMerged$D6d_n.log2FC > 1 & finalMerged$D6d_n.padj < 0.01) &
#                           (finalMerged$D1m_n.log2FC > 1 & finalMerged$D1m_n.padj < 0.01))
# )  &
#   (finalMerged$median.Development < finalMerged$median.Diapause)
# ,]
# length(dpUp1[,1]) # 2735
# write.table(dpUp1, file = "dpUp_strict.txt", row.names = F, sep = "\t")
# 
# 
# # List 1: Upregulated wrt both preD and nonD with FDR 0.05
# dpUp2 = finalMerged[(((finalMerged$D3d.log2FC > 0 & finalMerged$D3d.padj < 0.05) &
#                     (finalMerged$D6d.log2FC > 0 & finalMerged$D6d.padj < 0.05) &
#                     (finalMerged$D1m.log2FC > 0 & finalMerged$D1m.padj < 0.05)) |
#                     ((finalMerged$D3d_n.log2FC > 0 & finalMerged$D3d_n.padj < 0.05) &
#                     (finalMerged$D6d_n.log2FC > 0 & finalMerged$D6d_n.padj < 0.05) &
#                     (finalMerged$D1m_n.log2FC > 0 & finalMerged$D1m_n.padj < 0.05))
#                    )  &
#                      (finalMerged$median.Development < finalMerged$median.Diapause)
#                   ,]
# 
# length(dpUp2[,1]) # 5309
# write.table(dpUp2, file = "dpUp_intermediate.txt", row.names = F, sep = "\t")

# List 2: Upregulated wrt one of the preD and nonD with FDR 0.05 ********** FINAL USED FOR THE ANALYSIS **************
dpUp3 = finalMerged[( ( (finalMerged$D3d.log2FC > 0 & finalMerged$D3d.padj < 0.05) |
                        (finalMerged$D6d.log2FC > 0 & finalMerged$D6d.padj < 0.05) |
                        (finalMerged$D1m.log2FC > 0 & finalMerged$D1m.padj < 0.05)
                      ) |
                      ( (finalMerged$D3d_n.log2FC > 0 & finalMerged$D3d_n.padj < 0.05) |
                        (finalMerged$D6d_n.log2FC > 0 & finalMerged$D6d_n.padj < 0.05) |
                        (finalMerged$D1m_n.log2FC > 0 & finalMerged$D1m_n.padj < 0.05)
                      )
                    ) &
                      (finalMerged$median.Development < finalMerged$median.Diapause)
                    ,]

length(dpUp3[,1])
write.table(dpUp3, file = "2_Counts_DE/dpUp_relaxed.txt", row.names = F, sep = "\t")


################### DOWN ###############################
# List 2: Significantly down wrt both preD and nonD with FDR 0.05
# dpDown1 = finalMerged[(((finalMerged$D3d.log2FC < 0 & finalMerged$D3d.padj < 0.05) &
#                       (finalMerged$D6d.log2FC < 0 & finalMerged$D6d.padj < 0.05) &
#                       (finalMerged$D1m.log2FC < 0 & finalMerged$D1m.padj < 0.05)) |
#                       ((finalMerged$D3d_n.log2FC < 0 & finalMerged$D3d_n.padj < 0.05) &
#                       (finalMerged$D6d_n.log2FC < 0 & finalMerged$D6d_n.padj < 0.05) &
#                       (finalMerged$D1m_n.log2FC < 0 & finalMerged$D1m_n.padj < 0.05))
# )
# ,]
# length(dpDown1[,1]) # 5320
# write.table(dpDown1, file = "dpDown_strict.txt", row.names = F, sep = "\t")
# 
# # List 2: down wrt both preD and nonD
# dpDown2 = finalMerged[((finalMerged$D3d.log2FC < 0) &
#                          (finalMerged$D6d.log2FC < 0) &
#                          (finalMerged$D1m.log2FC < 0) &
#                          (finalMerged$D3d_n.log2FC < 0) &
#                          (finalMerged$D6d_n.log2FC < 0) &
#                          (finalMerged$D1m_n.log2FC < 0)
#  )
#  ,]
# 
# length(dpDown2[,1]) # 6992
# write.table(dpDown2, file = "dpDown_intermediate.txt", row.names = F, sep = "\t")

# Final used for paper ****************************************
dpDown3 = finalMerged[(finalMerged$median.Development > finalMerged$median.Diapause),]
length(dpDown3[,1])
write.table(dpDown3, file = "2_Counts_DE/dpDown_relaxed.txt", row.names = F, sep = "\t")

# Plot boxplots for most relaxed data
# boxplot(dpUp3$D3d.log2FC, dpUp3$D6d.log2FC, dpUp3$D1m.log2FC,dpDown3$D3d.log2FC, dpDown3$D6d.log2FC, dpDown3$D1m.log2FC, outline = F)
# boxplot(dpUp3$median.Diapause, dpUp3$median.Development, dpDown3$median.Diapause, dpDown3$median.Development, outline = F)
# 
# boxplot(dpUp3$median.Diapause, dpUp3$median.Development, outline = F)
# boxplot(dpDown3$median.Diapause, dpDown3$median.Development, outline = F)

# dpUp3$paralog = "up"
# dpDown3$paralog = "down"
# xData = rbind(dpUp3, dpDown3)
# 
# head(xData)
# 
# library(ggplot2)
# p1 <- ggplot(xData, aes(y = median.Diapause, x = paralog)) + 
#   geom_boxplot(outlier.shape = NA) +
#   ylim(0, 1000)
# p1

#p1 <- ggplot(xData, aes(y = median.Diapause, x = paralog)) + 
#  geom_violin(outlier.shape = NA) +
#  ylim(0, 1000)
#p1


# p1 + geom_jitter(shape=16, position=position_jitter(0.2))
# 
# p2 <- ggplot(xData, aes(y = median.Development, x = paralog)) + 
#   geom_boxplot(outlier.shape = NA,notch=T) +
#   ylim(0, 1000)
# p2

# # Get up and down based on median values
# dpUp = finalMerged[(finalMerged$median.FC.Diapause > 0),]
# length(dpUp[,1])
# write.table(dpUp, file = "dpUp_median_all.txt", row.names = F, sep = "\t")
# 
# dpUpS = finalMerged[(finalMerged$median.FC.Diapause > 0 & finalMerged$median.FC.padj < 0.05),]
# length(dpUpS[,1])
# write.table(dpUpS, file = "dpUp_median_significant.txt", row.names = F, sep = "\t")
# 
# 
# dpDown = finalMerged[(finalMerged$median.FC.Diapause < 0),]
# length(dpDown[,1])
# write.table(dpDown, file = "dpDown_median_all.txt", row.names = F, sep = "\t")
# 
# dpDownS = finalMerged[(finalMerged$median.FC.Diapause < 0 & finalMerged$median.FC.padj < 0.05),]
# length(dpDownS[,1])
# write.table(dpDownS, file = "dpDown_median_significant.txt", row.names = F, sep = "\t")
# 
# boxplot(dpUp$median.FC.Diapause, dpDown$median.FC.Diapause, outline = F)
# 
# hist(finalMerged$median.FC.Diapause, 1000)

