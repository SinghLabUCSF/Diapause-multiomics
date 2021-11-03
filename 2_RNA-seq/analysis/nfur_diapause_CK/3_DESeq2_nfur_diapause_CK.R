#  _ __  __ _ _ _ __ _ _ __  
# | '_ \/ _` | '_/ _` | '  \ 
# | .__/\__,_|_| \__,_|_|_|_|
# |_| 
# 
# To do differential enrichment analysis betweeen diapause and development

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

sampleTable <- read.csv("../../fastq/nfur_diapause_CK/experiment_design.csv", row.names = 1)  
sampleTable

library("DESeq2")

countdata = read.csv("1_Raw_counts/Counts_nfur_diapause_CK.csv", header = T, row.names = 1)
colnames(countdata)
coldata <- DataFrame(sampleTable)
rownames(coldata)
colnames(countdata) == rownames(coldata) # To make sure that the column names are identical

# Make dds object for DESEQ2
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

# Model with batch effects
#dds <- DESeqDataSetFromMatrix(countData = countdata,
#                              colData = coldata,
#                              design = ~ batch + condition)
  
# Check if condition is the correct column for ~
coldata

# only keep rows that have at least 15 sum count
dds <- dds[ rowSums(counts(dds)) > 15, ]

# Plot PCA using rld normalization - TAKES TIME ------------------------------------------------------
#rld <- rlog(dds, blind = FALSE)
#head(assay(rld), 3)
#plotPCA(rld, intgroup = c("condition"))

# Plot PCA using vst normalization ------------------------------------------------------
vsd <- vst(dds)
head(assay(vsd), 3)
plotPCA(vsd, intgroup = c("condition"))

dds <- DESeq(dds)

# Normalized counts
dds <- estimateSizeFactors(dds)
normcounts <- counts(dds, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable)
write.table(normcounts, file="2_Counts_DE/CountsNormDESeq2_nfur_diapause_CK.csv", sep= ",")
head(normcounts)

# Calculate differential expression. The reference (untretaed) is the second in contrast function for pairwise between each condition
# So for ("condition", "D3d","PreD"), PreD is untreated (reference) and D3d is treated (comparison)
# Fold changes will be wrt untreated samples, second in contrast
difD3d <- results (dds,contrast=c("condition", "D3d","PreD"))
write.table(difD3d, file="2_Counts_DE/DE_D3d-vs-preD.txt", sep="\t")

difD6d <- results (dds,contrast=c("condition", "D6d","PreD"))
write.table(difD6d, file="2_Counts_DE/DE_D6d-vs-preD.txt", sep="\t")

difD1m <- results (dds,contrast=c("condition", "D1m","PreD"))
write.table(difD1m, file="2_Counts_DE/DE_D1m-vs-preD.txt", sep="\t")

difNonD <- results (dds,contrast=c("condition", "NonD","PreD"))
write.table(difNonD, file="2_Counts_DE/DE_NonD-vs-preD.txt", sep="\t")

difD3d2 <- results (dds,contrast=c("condition", "D3d","NonD"))
write.table(difD3d2, file="2_Counts_DE/DE_D3d-vs-nonD.txt", sep="\t")

difD6d2 <- results (dds,contrast=c("condition", "D6d","NonD"))
write.table(difD6d2, file="2_Counts_DE/DE_D6d-vs-nonD.txt", sep="\t")

difD1m2 <- results (dds,contrast=c("condition", "D1m","NonD"))
write.table(difD1m2, file="2_Counts_DE/DE_D1m-vs-nonD.txt", sep="\t")

# # Read length file and make TPM from normalized counts - I noticed that the TPMs are almost identical before and after normalization though
# len = read.csv("nfur_diapause_CK/GeneLengths_nfur_diapause_CK.csv", header = T, row.names = 1)
# head(normcounts)
# fortpm = merge(normcounts,len[1], by="row.names")
# fortpm <- data.frame(koko[,-1], row.names=koko[,1])
# head(fortpm)
# 
# # using sapply to apply the tpm function on the data
# normtpm = sapply(fortpm[1:15], tpm, length = fortpm$PreD1.y)
# row.names(normtpm) = row.names(fortpm)
# head(normtpm)
# apply(normtpm, 2, sum) # sums should be 1 million for each column
# write.table(normtpm, file="nfur_diapause_CK/TPMNormalized_nfur_diapause_CK.csv", sep= ",", row.names = T)
