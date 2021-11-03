# Generate a combined plot for Figure 2 for Alim and Nfur comparison
library("ggpubr") # To combine panels
library(ggplot2)
library(eulerr) # for venn

# Set wd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file = read.table("nfurzeri_alimnaeus_paralogs.txt", header = T, sep = "\t")
file = file[complete.cases(file), ]
head(file)

pdf(file ="Alim_nfur_comparison.pdf", height = 5, width = 5)
par(mfrow=c(1,2))
boxplot(file$median.Diapause.1[file$nfur.category == 'NeoF'],
        file$median.Development.1[file$nfur.category == 'NeoF'],
        file$median.Diapause.alim.1[file$nfur.category == 'NeoF'],
        file$median.Development.alim.1[file$nfur.category == 'NeoF'],
        outline = F)

boxplot(file$median.Diapause.2[file$nfur.category == 'NeoF'],
        file$median.Development.2[file$nfur.category == 'NeoF'],
        file$median.Diapause.alim.2[file$nfur.category == 'NeoF'],
        file$median.Development.alim.2[file$nfur.category == 'NeoF'],
        outline = F)
dev.off()

# boxplot(file$median.Diapause.1[file$nfur.category == 'Unclassified'],
#         file$median.Development.1[file$nfur.category == 'Unclassified'],
#         file$median.Diapause.alim.1[file$nfur.category == 'Unclassified'],
#         file$median.Development.alim.1[file$nfur.category == 'Unclassified'],
#         outline = F)

#file = file[file$nfur.category == 'NeoF',]
#ggplot(file, aes(x=file$median.FC.Diapause.1, y=file$Dia.log2FC.alim.1)) + geom_point() +geom_smooth(method=lm)
#ggplot(file, aes(x=file$median.Diapause.1, y=file$median.Diapause.alim.1[file$nfur.category == 'NeoF'])) + geom_point() +geom_smooth(method=lm)

# Some correlations that I am not plotting right now
# ggplot(file, aes(x=log(file$median.Diapause.1), y=log(file$median.Diapause.alim.1))) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Diapause.2), y=log(file$median.Diapause.alim.2))) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Development.1), y=log(file$median.Development.alim.1))) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Development.2), y=log(file$median.Development.alim.2))) + geom_point() +geom_smooth(method=lm)
# 
# ggplot(file, aes(x=log(file$median.Diapause.1), y=log(file$median.Diapause.alim.2))) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Diapause.2), y=log(file$median.Diapause.alim.1))) + geom_point() +geom_smooth(method=lm)
 
#ggplot(file, aes(x=file$median.Diapause.1, y=file$median.Diapause.alim.1[file$nfur.category == 'NeoF'])) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Diapause.1), y=log(file$median.Development.alim.1))) + geom_point() +geom_smooth(method=lm)
# ggplot(file, aes(x=log(file$median.Diapause.1), y=log(file$median.Diapause.2))) + geom_point() +geom_smooth(method=lm)

# Generate correlation. First get all genes with at least 10 tpm for all samples, and significant in both ALim and Nfur with FDR 0.01
file2 = file[(file$median.FC.padj.1 < 0.01 & file$Dia.padj.alim.2 < 0.01 &
                      file$median.Development.1 >= 10 &
                      file$median.Development.2 >= 10 &
                      file$median.Diapause.1 >= 10 &
                      file$median.Diapause.2 >= 10 &
                      file$median.Diapause.alim.1 >= 10 &
                      file$median.Development.alim.1 >= 10 &
                      file$median.Development.alim.2 >= 10 &
                      file$median.Diapause.alim.2 >= 10 
                      ),]

# ggplot(file2, aes(x=file2$median.FC.Diapause.1,
#                  y=file2$Dia.log2FC.alim.1)) +
#         geom_point() +
#         geom_smooth(method=lm) 

# Alim Nfur fold change correlation
pdf(file ="Alim_nfur_FC_correlation.pdf", height = 5, width = 4)

plot (file2$median.FC.Diapause.1, file2$Dia.log2FC.alim.1, 
      main = paste ("FC: Spearman ", cor(file2$median.FC.Diapause.1, file2$Dia.log2FC.alim.1, method = "spearman"), "; p ", cor.test (file2$median.FC.Diapause.1, file2$Dia.log2FC.alim.1)$p.value),
      pch = 16,
      col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3),
      xlab = "N. furzeri", ylab = "A. limnaeus"
)
abline(lm(file2$median.FC.Diapause.1 ~ file2$Dia.log2FC.alim.1))
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

dev.off()

# Overlap in NeoF genes
# test for enrichment of NeoF in Alim and Nfur 
phyper((606-1), 3417, (9854-3417), 1341, lower.tail = F, log.p = FALSE)
# 6.787643e-18

pdf(file ="Alim_nfur_NeoFOverlap.pdf", height = 5, width = 4)
NeoF_venn <-  c(Nfur = 3417, Alim = 1341,
                "Nfur&Alim" = 606)
fit3 <- euler(NeoF_venn, shape = "ellipse")
plot(fit3, quantities = TRUE)

dev.off()




