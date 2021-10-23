# THIS WAS RUN FOR PAPER USING DIFFBIND VERSION 2.16.2
library("DiffBind")
library("tidyverse")

# Set wd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############# INPUT FILES #################
# Nfur
sample_sheet = 'nfur_EdgeR/sample_sheet_3.csv'
output_dir = 'nfur_EdgeR/'
# Medaka
#sample_sheet = 'medaka/sample_sheet_medaka.csv'
#output_dir = 'medaka/'
# Medaka - selected stages
#sample_sheet = 'olat_selected-stages/sample_sheet_medaka.csv'
#output_dir = 'olat_selected-stages/'
# Alim
#sample_sheet = 'alim/sample_sheet_alim.csv'
#output_dir = 'alim/'
# Drer
#sample_sheet = 'drer/sample_sheet_drer.csv'
#output_dir = 'drer/'
# Drer selected stages
#sample_sheet = 'drer_selected-stages/sample_sheet_drer.csv'
#output_dir = 'drer_selected-stages/'

# All aligned to nfur
#sample_sheet = 'all_aligned_to_nfur/sample_sheet_killi.csv'
#output_dir = 'all_aligned_to_nfur/'
###########################################

# Read sample description sheet. See help for the dba function for column specification
# Column names are fixed so it automatically reads the peaks etc.
samples <- read.csv(sample_sheet) 
samples

# Load and read the sample sheet, peaks and metadata
DBdata <- dba(sampleSheet=samples)
DBdata

# plot a correlation heatmap that gives an initial clustering
# of the samples using cross-correaltions of the rows in the 
# binding matrix.
#plot(DBdata) 

# calculate a binding matrix based on the read counts for every sample
# WARNING: This also reads bam files so it will take a while to run. If you just want to check quickly, change csv file to just have one condition.
DBdata <- dba.count(DBdata) 
DBdata

# group samples based on condition for contrast 
# Also have to set the minMembers parameter to 2 
DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, minMembers = 2)
DBdata

# run the main function (dba.analyze) of DiffBind that performs differential binding
# Use both methods for comparison
DBdata <- dba.analyze(DBdata, method = DBA_DESEQ2)
DBdata <- dba.analyze(DBdata, method = DBA_EDGER)
DBdata

# **** Extract the consensus peak set for all samples and the RPKM counts for each condition
# Remember I am not subtracting 0 from start so BED format specifications may not meet strictly. This is to be consistent in all the DESEq2 scripts. Rething and revisit later.
# I rearranaged columns in the beginning to make it competible to bed format
# See here https://www.biostars.org/p/269615/
consensusPeaks = dba.peakset(DBdata, peak.format = "bed", consensus = T, bRetrieve = T)
consensusPeaksDF = data.frame(chromosome=seqnames(consensusPeaks),
                              starts=start(consensusPeaks), # BED uses 0-based coordinates so I can do -1 but for now I am keeping this
                              ends=end(consensusPeaks),
                              names=paste0("peak_", names(consensusPeaks)),
                              scores=c(rep("0", length(consensusPeaks))),
                              strand=strand(consensusPeaks),
                              consensusPeaks)
write.table(consensusPeaksDF, file = paste0(output_dir, "consensusPeaks_rpkm.txt"), sep = "\t", quote = F, row.names = F) 


# Extract the RPKM, read counts etc. for each sample
# This is same as above but has more information.
for (i in as.numeric(rownames(dba.show(DBdata)))){ # I am getting the rownames for each sample that runs from 1 to 10
  
  print (i)
  print (DBdata$samples$SampleID[i]) # get the sample id for the name of file
  sampleid = DBdata$samples$SampleID[i]
  countfilename = paste0(output_dir, sampleid, "_Readcounts.txt")
  
  # DBdata$peaks[1:10] have the count matrices that I need, so print them to files called read counts
  write.table(DBdata$peaks[i], file = countfilename, sep = "\t", quote = F, row.names = F)
}
rm(i)

# Output the number of peaks 
DEcounts = dba.show(DBdata, bContrasts = T)
write.csv(DEcounts, file = paste0(output_dir, "DEcounts.csv"), quote = F)

# Plot heatmap and bindinf affinity correation using only the differentially bound sites.
pdf(paste0(output_dir, "Correlation_matrix.pdf"))
dba.plotHeatmap(DBdata, correlations=T, method = DBA_EDGER)
dev.off()

pdf(paste0(output_dir, "Binding-affinity_EdgeR.pdf"))
dba.plotHeatmap(DBdata, correlations=F, scale = "row", method = DBA_EDGER)
dev.off() 

# Global PCA
pdf(paste0(output_dir, "Global_PCA.pdf"))
dba.plotPCA(DBdata)
dev.off()

# Subset the DBdata object for plotting etc
#msk = dba.mask(DBdata, "Condition", c("Dia6D", "Escape", "Dia1m", "Kvyoung", "Alim_1M_1", "Alim_Dev_1", "Alim_KV_1"), combine='or')
msk = dba.mask(DBdata, "Condition", c("Escape", "Dia1m", "Alim_1M_1", "Alim_Dev_1"), combine='or')
db2  = dba(DBdata, mask=msk)
dba.plotPCA(db2)

# Do specific stuff for contrasts -- using ***EDGER*** as method here
for (i in 1:length(DEcounts[,1])){
  
  filename = paste0(output_dir, DEcounts[i,1], "-", DEcounts[i,3]) 
  print (filename)
  
  # Plot and save contrast specific files 
  # Specific PCA
  pdf(paste0(filename, "_PCA_EdgeR.pdf"))
  dba.plotPCA(DBdata, contrast = i, method = DBA_EDGER)
  dev.off()

  # Boxplots for the read distribution differences between the classes for selected contrasts.
  # The first two boxes show distribution of reads over all differentially bound sites; 
  # the middle two show differences on those sites where the affinity increases and the two boxes on the right show differences where the affinity increases in
  # highlighted samples 
  pdf(paste0(filename, "_Read-distribution_EdgeR.pdf"))
  dba.plotBox(DBdata, contrast = i, method = DBA_EDGER)
  dev.off() 
  
  # display the binding affinity heatmap to see the binding patterns across the identified regions. You can control
  pdf(paste0(filename, "_Binding-affinity_EdgeR.pdf"))
  dba.plotHeatmap(DBdata, contrast=i, correlations=F, method = DBA_EDGER)
  dev.off() 

  pdf(paste0(filename, "_Correlation_EdgeR.pdf"))
  dba.plotHeatmap(DBdata, contrast=i, correlations=T, method = DBA_EDGER)
  dev.off() 
  
  # MA plots
  pdf(paste0(filename, "_MA_EdgeR.pdf"))
  dba.plotMA(DBdata, contrast = i, method = DBA_EDGER)
  dev.off()
  
  # Volcano plots
  pdf(paste0(filename, "_Volcano_EdgeR.pdf"))
  dba.plotVolcano(DBdata, contrast = i, method = DBA_EDGER)  
  dev.off()
  
  # *********************************************
  # report the differentially bound peak regions, identified by either method (DESeq2/edgeR)
  report = dba.report(DBdata, contrast = i, method=DBA_EDGER)

  # Export 3 bed files - up, down and all DE sites
  # First make a data frame with column order as in bed file. Based on https://www.biostars.org/p/89341/
  df1 <- data.frame(seqnames=seqnames(report),
                    starts=start(report)-1, # BED uses 0-based coordinates so do -1 from start
                    ends=end(report),
                    names=paste0("peak_", names(report)), # Peak names. These are consistent for all consensus peaks
                    scores=c(rep("0", length(report))),
                    strand=c(rep(".", length(report))),
                    fold=elementMetadata(report)$Fold
  )
  # Then write 3 bed files. Select columns 1:6 and avoid fold change. you can uncomment for checking 
  write.table(df1[df1$fold > 0,1:6], file = paste0(filename, "_DE-up_EdgeR.bed"), sep = "\t", quote = F, row.names = F)
  write.table(df1[df1$fold < 0,1:6], file = paste0(filename, "_DE-down_EdgeR.bed"), sep = "\t", quote = F, row.names = F)
  write.table(df1[,1:6], file = paste0(filename, "_DE-all_EdgeR.bed"), sep = "\t", quote = F, row.names = F)
  
  # Make another data frame to get additional columns, FDR, p-values etc.
  df2 = data.frame(report)
  df2$start = df2$start - 1 # BED uses 0-based coordinates so do -1 from start
  
  # Write a csv file with all columns  
  finaldf = data.frame(df1[,1:6], df2[,c(4,6:11)])
  write.table(finaldf, file = paste0(filename, "_DE-all_EdgeR.csv"), sep = ",", quote = F, row.names = F)
  
  rm(filename, report, df1, df2, finaldf) # remove variables for next run
}


