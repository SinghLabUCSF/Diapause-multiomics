# Setwd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

dia = read.table("RNA-ATAC_correlation_diapause_all.txt", sep = "\t", header = T)
dev = read.table("RNA-ATAC_correlation_development_all.txt", sep = "\t", header = T)

head(dia)
head(dev)

# Remove RNA-seq. ATAC-seq is never 0 in this data
dia = dia[dia$rna > 0,]
dev = dev[dev$rna > 0,]

# DIAPAUSE - MEAN ATAC correlation
# pdf(file = "RNA-ATAC-correlation_diapause.pdf")
# ggplot(dia, aes(x=log(rna), y=log(atac_median))) + 
#   geom_point(alpha = 0.3) + 
#   geom_smooth(method=lm) + 
#   theme_bw()
# dev.off()
# cor.test(dia$rna, dia$atac_median, method = "spearman")

pdf(file = "RNA-ATAC-correlation_diapause.pdf")
plot (log(dia$rna) ~ log(dia$atac_median), 
      main = paste ("FC: Spearman ", cor(log(dia$rna), log(dia$atac_median), method = "spearman"), "; p ", cor.test (log(dia$rna), log(dia$atac_median))$p.value),
      pch = 16,
      col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3),
      xlab = "log (median ATAC-seq accessibility)", ylab = "log(median RNA-seq expression)"
)
abline(lm(log(dia$rna) ~ log(dia$atac_median)), untf = T)
dev.off()


# DEVELOPMENT - MEAN ATAC
head(dev)
# pdf(file = "RNA-ATAC-correlation_development.pdf")
# ggplot(dev, aes(x=log(rna), y=log(atac_median))) + 
#   geom_point(alpha = 0.3) + 
#   geom_smooth(method=lm) +
#   theme_bw()
# dev.off()
# 
# cor.test(dev$rna, dev$atac_median, method = "spearman")

# This may have reduced file size but line is not plotting properly
pdf(file = "RNA-ATAC-correlation_development.pdf")
plot (log(dev$rna) ~ log(dev$atac_median), 
      main = paste ("FC: Spearman ", cor(log(dev$rna), log(dev$atac_median), method = "spearman"), "; p ", cor.test (log(dev$rna), log(dev$atac_median))$p.value),
      pch = 16,
      col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3),
      xlab = "log (median ATAC-seq accessibility)", ylab = "log(median RNA-seq expression)"
)
abline(lm(log(dev$rna) ~ log(dev$atac_median)), untf = T)
dev.off()

#rm(dia,dev)

# # Correlation with only promoter peak. This is also significant but I am not using it
# diaP = read.table("RNA-ATAC_correlation_diapause_promoter.txt", sep = "\t", header = T)
# devP = read.table("RNA-ATAC_correlation_development_promoter.txt", sep = "\t", header = T)
# 
# diaP = diaP[diaP$rna > 0,]
# devP = devP[devP$rna > 0,]
# 
# ggplot(diaP, aes(x=log(rna), y=log(atac_median))) + geom_point() + geom_smooth(method=lm)
# cor.test(diaP$rna, diaP$atac_median, method = "spearman")
# 
# ggplot(devP, aes(x=log(rna), y=log(atac_median))) + geom_point() + geom_smooth(method=lm)
# cor.test(devP$rna, devP$atac_median, method = "spearman")
# 
