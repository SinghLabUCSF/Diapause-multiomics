# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(UpSetR)
library(ggplot2)

# Intersections of all sets 
# peaks = read.table(file = "/Volumes/Mybook_3/ATAC_Seq_Killifish/Analysis/nfur/MasterPeakFile_UP-all.txt", sep = "\t", header = T, row.names = 1)
# head(peaks)
# upset(peaks, nsets = 10, nintersects = 10000) # Total intersections are  lot here
# 
# # Intersection of all diapause sets
# peaks = read.table(file = "/Volumes/Mybook_3/ATAC_Seq_Killifish/Analysis/nfur/MasterPeakFile_UP-diapause.txt", sep = "\t", header = T, row.names = 1)
# head(peaks)
# upset(peaks, nsets = 6, nintersects = 100)

# Intersection of all diapause and development sets
peaks = read.table(file = "MasterPeakFile_DEseq2-EdgeR_Combined_OnlyDE_UP.txt", sep = "\t", header = T, row.names = 1)
dim(peaks)
head(peaks)
upset(peaks, nsets = 5, nintersects = 100)
