# To do differential enrichment analysis by combining multiple samples

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

sampleTable <- read.csv("/Volumes/Mybook_3/RNA-seq/fastq/alim/dia_nondia/experiment_design_D-nonD.csv", row.names = 1)
sampleTable

library("DESeq2")

countdata = read.csv("Counts_alim_dia-nondia.csv", row.names = 1)
colnames(countdata)
coldata <- DataFrame(sampleTable)
rownames(coldata)
colnames(countdata) == rownames(coldata) # To make sure that the column names are identical

# Make dds object for DESEQ2
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

# Model with batch effects
# dds <- DESeqDataSetFromMatrix(countData = countdata,
#                               colData = coldata,
#                               design = ~ batch + condition)
  
# Check if condition is the correct column for ~
coldata

# only keep rows that have at least 15 sum count
dds <- dds[ rowSums(counts(dds)) > 15, ]

# Plot PCA using rld normalization (takes longer so I will use vsd) ------------------------------------------------------
# rld <- rlog(dds, blind = FALSE)
# head(assay(rld), 3)
# plotPCA(rld, intgroup = c("condition"))

# Plot PCA using vst normalization ------------------------------------------------------
vsd <- vst(dds)
head(assay(vsd), 3)
z = plotPCA(vsd, intgroup = c("condition"))
z
# to edit PCA a bit
library(ggplot2)
library(ggrepel)
colnames(sampleTable)
z + geom_text_repel(aes(label = sampleTable$condition))
z + geom_text_repel(aes(label = rownames(sampleTable)))

dds <- DESeq(dds)

# Normalized counts
dds <- estimateSizeFactors(dds)
normcounts <- counts(dds, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable)
write.table(normcounts, file="CountsNormDESeq2_alim_dia-nondia.csv", sep= ",")
head(normcounts)

# Calculate differential expression. The reference (untretaed) is the second in contrast function
# So for ("condition", "D3d","PreD"), PreD is untreated (reference) and D3d is treated (comparison)
# Fold changes will be wrt untreated samples i.e. second in contrast
difDia <- results (dds,contrast=c("condition", "AlimDia","Alim4dpD"))
write.table(difDia, file="DE_Dia-vs-4dpD.txt", sep="\t")


