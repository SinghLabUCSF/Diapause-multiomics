# To do differential enrichment analysis by combining multiple samples

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("DESeq2")

## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/
# rpkm <- function(counts, lengths) {
#   rate <- counts / lengths
#   rate / sum(counts) * 1e9
# }
# 
# tpm <- function(counts, lengths) {
#   rate <- counts / lengths
#   rate / sum(rate) * 1e6
# }

sampleTable <- read.csv("../../../fastq/nfur_ast_aaul/experiment_design_nfur_aaul_ast_combined.csv", row.names = 1)
sampleTable

countdata = read.csv("nfur_aaul-ast_ortholog_count_matrix.csv", row.names = 1)
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
# z + geom_text_repel(aes(label = sampleTable$condition))
# z + geom_text_repel(aes(label = sampleTable$original_names))
# z + geom_text_repel(aes(label = sampleTable$prefix))
z + geom_text_repel(aes(label = sampleTable$lib))

# PCA for selected subgroups
# vsd$condition
# vsd.sub <- vsd[ , vsd$condition %in% c("nfur_kv_young", "aaul_kv", "ast_kv_young")]

# pp = plotPCA_New(vsd.sub, intgroup=c("condition"), ax1="PC1", ay1="PC2", av1 = 1, av2 = 2)
# pp + geom_text_repel(aes(label = vsd.sub$condition))


dds <- DESeq(dds)

# Normalized counts
dds <- estimateSizeFactors(dds)
normcounts <- counts(dds, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable)
write.table(normcounts, file="CountsNormDESeq2_aaul-ast-nfur.csv", sep= ",")
head(normcounts)

# # Calculate differential expression. The reference (untretaed) is the second in contrast function
# # So for ("condition", "D3d","PreD"), PreD is untreated (reference) and D3d is treated (comparison)
# # Fold changes will be wrt untreated samples i.e. second in contrast
# difD3d <- results (dds,contrast=c("condition", "nfur_kv_young","nfur_kv_old"))
# write.table(difD3d, file="DE_D3d-vs-preD.txt", sep="\t")
# 
# difD6d <- results (dds,contrast=c("condition", "D6d","PreD"))
# write.table(difD6d, file="DE_D6d-vs-preD.txt", sep="\t")
# 
# difD1m <- results (dds,contrast=c("condition", "D1m","PreD"))
# write.table(difD1m, file="DE_D1m-vs-preD.txt", sep="\t")
# 
# difNonD <- results (dds,contrast=c("condition", "NonD","PreD"))
# write.table(difNonD, file="DE_NonD-vs-preD.txt", sep="\t")
# 
# difD3d2 <- results (dds,contrast=c("condition", "D3d","NonD"))
# write.table(difD3d2, file="DE_D3d-vs-nonD.txt", sep="\t")
# 
# difD6d2 <- results (dds,contrast=c("condition", "D6d","NonD"))
# write.table(difD6d2, file="DE_D6d-vs-nonD.txt", sep="\t")
# 
# difD1m2 <- results (dds,contrast=c("condition", "D1m","NonD"))
# write.table(difD1m2, file="DE_D1m-vs-nonD.txt", sep="\t")

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




plotPCA_New <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, ax1 = "PC1", ay1 = "PC2", av1 = "1", av2 = "2") 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = ax1, y = ay1, color = "group")) + 
    geom_point(size = 3) + xlab(paste0(ax1, ": ", round(percentVar[av1] * 
                                                          100), "% variance")) + ylab(paste0(ay1, ": ", round(percentVar[av2] * 100), "% variance")) + coord_fixed()
}

