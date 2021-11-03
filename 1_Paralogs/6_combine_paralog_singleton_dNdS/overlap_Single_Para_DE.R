# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Get all combined paralogs
paraFile = read.table("All_paraogs_20200820.txt", head = T, sep = "\t")
head(paraFile)
# Get unambiguous paralogs - that are in all the comparisons
paraFile = paraFile[(paraFile$X71spp != '' & paraFile$X31spp != '' & paraFile$X13spp != '' & paraFile$Ensembl != ''),]
paralogs = unique(c(paraFile$Id1, paraFile$Id2)) # get unique gene ids
length(paralogs)

# Get all combined singletons
singleFile = read.table("All_singletons_20210211.txt", head = T, sep = "\t")
head(singleFile)
# get unambiguous singletons that are in all the comparisons
singleFile = singleFile[(singleFile$X71spp == 'Y' & singleFile$X31spp == 'Y' & singleFile$X13spp == 'Y' & singleFile$Ensembl == 'Y'),]
singletons = singleFile$Id # get unique singleton gene ids
length (singletons)

# Remove singletons and paralogs that overlap - ideally none
# Final paralogs
finalPara = paralogs[!(paralogs %in% singletons)]
paste0("Total paralogs: ", length(finalPara))
totalPara = length(finalPara)
print (totalPara)

# Final singletons
finalSingle = singletons[!(singletons %in% paralogs)]
paste0("Total singletons: ", length(finalSingle))
totalSingle = length(finalSingle)
print (totalSingle)

# Get all DE up in diapause
DEfile = read.table("../../RNA-seq/analysis/nfur_diapause_CK/2_Counts_DE/dpUp_relaxed.txt", head = T)
head(DEfile)
differential = DEfile$genes
head(differential)

# Final DE - remove the ones not in any
finalDE = differential[(differential %in% singletons) | (differential %in% paralogs)]
paste0("Total differential genes: ", length(finalDE))

paste0("Total DE genes in singletons: ", length(finalDE [finalDE  %in% finalSingle]))
singleDE = length(finalDE [finalDE  %in% finalSingle])
paste0("Total DE genes in paralogs: ", length(finalDE [finalDE  %in% finalPara]))
paraDE = length(finalDE [finalDE  %in% finalPara])

# Chi-square test for the distribution of DE in paralog and singletons
chisq.test(c(paraDE, singleDE), p = c(totalPara/(totalPara+totalSingle), totalSingle/(totalPara+totalSingle)))

# Bar chart of DE genes in paralogs and singletons
pdf(file = "Singleton_Paralog_distribution.pdf")
par(mar=c(8,5,3,3))
barplot(
  c(totalPara*100/(totalPara+totalSingle),
    paraDE*100/(paraDE+singleDE), 
    0,
    totalSingle*100/(totalPara+totalSingle),
    singleDE*100/(paraDE+singleDE)
  ), names.arg = c("Total Paralogs", "DE paralogs", "", "Total singletons", "DE singletons"), las = 2,
  ylab="Percentage of genes"
)
dev.off()

# # Pie chart of DE genes in paralogs and singletons - NOT USED
# slices <- c(length(finalDE [finalDE  %in% finalSingle]),
#             length(finalDE [finalDE  %in% finalPara]),
#             length(finalDE[!((finalDE  %in% finalPara) | (finalDE  %in% finalSingle))]))
# lbls <- c(paste0 ("DE with no paralogs: ", round(100*length(finalDE [finalDE  %in% finalSingle])/length(finalDE), 1), "%"),
#           paste0 ("DE with paralogs: ", round(100*length(finalDE [finalDE  %in% finalPara])/length(finalDE), 1), "%")
#           )
# pie(slices, labels = lbls)
