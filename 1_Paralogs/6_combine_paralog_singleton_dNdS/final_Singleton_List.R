# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Get all combined singletons
singleFile = read.table("All_singletons_20210211.txt", head = T, sep = "\t")
head(singleFile)

# Final singletons
singletons = singleFile[(singleFile$X71spp == 'Y' & singleFile$X31spp == 'Y' & singleFile$X13spp == 'Y' & singleFile$Ensembl == 'Y'),]
#singletons = singleFile[(singleFile$X71spp == 'Y'),]
finalSingle = singletons$Id
head(finalSingle)
paste0("Total singletons: ", length(finalSingle))

# Get all DE up in diapause
DEfile = read.table("../../RNA-seq/analysis/nfur_diapause_CK/2_Counts_DE/dpUp_relaxed.txt", head = T)
head(DEfile)
differential = DEfile$genes
head(differential)

# get single in Up diapause vs not 
singleUp = finalSingle[(finalSingle %in% differential)]
paste0("Total differential singleton genes: ", length(singleUp))

singleNoUp = finalSingle[!(finalSingle %in% differential)]
paste0("Total NOT differential singleton genes: ", length(singleNoUp))

write.table(finalSingle, "Singletons_final_All.txt", quote = F, row.names = F, col.names = F)
write.table(singleUp, "Singletons_final_UpInDiapuse.txt", quote = F, row.names = F, col.names = F)
write.table(singleNoUp, "Singletons_final_NotUpInDiapause.txt", quote = F, row.names = F, col.names = F)
