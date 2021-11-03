library("ggpubr") # To combine panels

# Set wd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("reshape2")
library("ggplot2")

file = read.table("4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt", header = T, sep = "\t")
head(file)

# I am doing KS test one sided, because wth 2 sides the exact p-value can not be calculated
p1 = ks.test(file$median.Diapause.1[file$Category == 'NeoF'], file$median.Development.1[file$Category == 'NeoF'], alternative = "l")$p.value
p2 = ks.test(file$median.Diapause.2[file$Category == 'NeoF'], file$median.Development.2[file$Category == 'NeoF'], alternative = "g")$p.value

pdf(file = "BoxPlot_1C.pdf", width = 6, height = 6)
boxplot(file$median.Diapause.1[file$Category == 'NeoF'],
        file$median.Development.1[file$Category == 'NeoF'],
        file$median.Diapause.2[file$Category == 'NeoF'],
        file$median.Development.2[file$Category == 'NeoF'],
        outline = F)
text(x = 1.5, y = 4000, labels = p1)
text(x = 3.5, y = 4000, labels = p2)
dev.off()

# boxplot(file$median.FC.Diapause.1[file$Category == 'NeoF'],
#         file$median.FC.Diapause.2[file$Category == 'NeoF'],
#         outline = F)
# 
# boxplot(file$median.Diapause.1[file$Category == 'NeoF'],
#          file$median.Development.1[file$Category == 'NeoF'],
#          outline = F)
#  
# boxplot(file$median.Diapause.2[file$Category == 'NeoF'],
#          file$median.Development.2[file$Category == 'NeoF'],
#          outline = F)
# 
# boxplot(file$median.Diapause.1[file$Category == 'NeoF'],
#         file$median.Diapause.2[file$Category == 'NeoF'],
#         outline = F)
# 
# boxplot(file$median.Development.1[file$Category == 'NeoF'],
#         file$median.Development.2[file$Category == 'NeoF'],
#         outline = F)
# 
# boxplot(file$median.Development.1[file$Category == 'Both down'],
#         file$median.Development.2[file$Category == 'Both down'],
#         outline = F)
# 
# boxplot(file$median.Development.1[file$Category == 'Unclassified'],
#         file$median.Development.2[file$Category == 'Unclassified'],
#         outline = F)
# 
# 
# boxplot(log((file$median.Diapause.1[file$Category == 'NeoF']/file$median.Diapause.2[file$Category == 'NeoF'])),
#         log((file$median.Diapause.1[file$Category == 'Unclassified']/file$median.Diapause.2[file$Category == 'Unclassified'])),
#         outline = F)

# boxplot(file$median.Diapause.1[file$Category == 'Both up'],
#         file$median.Development.1[file$Category == 'Both up'],
#         file$median.Diapause.2[file$Category == 'Both up'],
#         file$median.Development.2[file$Category == 'Both up'],
#         outline = F)
# 
# boxplot(file$median.Diapause.1[file$Category == 'Both down'],
#         file$median.Development.1[file$Category == 'Both down'],
#         file$median.Diapause.2[file$Category == 'Both down'],
#         file$median.Development.2[file$Category == 'Both down'],
#         outline = F)

# boxplot(file$median.FC.Diapause.1[file$Category == 'Both up'],
#         file$median.FC.Diapause.2[file$Category == 'Both up'],
#         outline = F)
# 
# boxplot(file$median.FC.Diapause.1[file$Category == 'Both down'],
#         file$median.FC.Diapause.2[file$Category == 'Both down'],
#         outline = F)
# 
# boxplot(file$median.FC.Diapause.1[file$Category == 'Unclassified'],
#         file$median.FC.Diapause.2[file$Category == 'Unclassified'],
#         outline = F)
# 
# 
# boxplot(file$median.FC.Diapause.1[file$Category == 'NeoF'],
#         file$median.FC.Diapause.2[file$Category == 'NeoF'],
#         file$median.FC.Diapause.1[file$Category == 'Both up'],
#         file$median.FC.Diapause.2[file$Category == 'Both up'],
#         file$median.FC.Diapause.1[file$Category == 'Both down'],
#         file$median.FC.Diapause.2[file$Category == 'Both down'],
#         file$median.FC.Diapause.1[file$Category == 'Unclassified'],
#         file$median.FC.Diapause.2[file$Category == 'Unclassified'],
#         outline = F)
# 
# colnames(file)
# #file2 = file[, c("gene1", "gene2" , "duplication_node", "Category", "median.FC.Diapause.1","median.FC.Diapause.2")]
# #colnames(file2)
# file2 = file[, c("gene1", "gene2" , "duplication_node", "Category", "median.Development.1", "median.Diapause.1", "median.Development.2", "median.Diapause.2")]
# #file4 = file[, c("gene1", "gene2" , "duplication_node", "Category", "median.Development.1", "median.Diapause.1")]
# #file5 = file[, c("gene1", "gene2" , "duplication_node", "Category", "median.Development.2", "median.Diapause.2")]
# colnames(file2)
# 
# file.m <- melt(file2)
# dim(file)
# head(file.m)
# 
# # file.m3 <- melt(file3)
# # file.m4 <- melt(file4)
# # file.m4 <- melt(file4)
# # head(file.m3)
# 
# ggplot(file.m, aes(x=factor(Category), y=value, fill = variable)) + 
#         #geom_boxplot(outlier.shape = NA,notch=T) 
#         geom_violin(trim=TRUE) 
#         #geom_boxplot(outlier.shape = NA, width=0.05) 
#         #geom_jitter(color="gray", size=0.4, alpha=0.9)
#         #geom_boxplot(outlier.shape = NA) 
# 
# ggplot(file.m[file.m$Category == "neoF",], aes(x=factor(Category), y=value, fill = variable)) + 
#         geom_violin() +
#         geom_boxplot(outlier.shape = NA, width=0.05) +
#         theme_classic()
# 
# `%nin%` <- Negate(`%in%`)
# 
# filteredList = file.m[file.m$Category == "NeoF",]
# head(filteredList)
# outliers1 = boxplot(filteredList[filteredList$variable == "median.Diapause.1", 6], plot = FALSE)$out
# outliers2 = boxplot(filteredList[filteredList$variable == "median.Development.1", 6], plot = FALSE)$out
# outliers3 = boxplot(filteredList[filteredList$variable == "median.Diapause.2", 6], plot = FALSE)$out
# outliers4 = boxplot(filteredList[filteredList$variable == "median.Development.2", 6], plot = FALSE)$out
# 
# head(outliers1)
# head(outliers2)
# head(outliers3)
# head(outliers4)
# 
# filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2 & filteredList$value %nin% outliers3 & filteredList$value %nin% outliers4),]
# 
# a = ggplot(filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2 & filteredList$value %nin% outliers3 & filteredList$value %nin% outliers4),], aes(x=factor(Category), y=value, fill = variable)) + 
#         geom_boxplot(outlier.shape = NA, width=0.05) +
#         geom_violin() +
#         geom_hline(yintercept = 0, colour="gray", linetype="dashed") +
#         theme_classic() 
# a
# 
# 
# rm(list = c("filteredList", "outliers1", "outliers2"))
# 
# filteredList = file.m[file.m$Category == "Both up",]
# head(filteredList)
# outliers1 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.1", 6], plot = FALSE)$out
# outliers2 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.2", 6], plot = FALSE)$out
# head(outliers1)
# head(outliers2)
# filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),]
# 
# rm(list = c("filteredList", "outliers1", "outliers2"))
# 
# 
# b = ggplot(filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),], aes(x=factor(Category), y=value, fill = variable)) + 
#         geom_boxplot(outlier.shape = NA, width=0.05) +
#         geom_violin() +
#         geom_hline(yintercept = 0, colour="gray", linetype="dashed") +
#         ylim(-3, +3) +
#         theme_classic() 
# rm(list = c("filteredList", "outliers1", "outliers2"))
# 
# filteredList = file.m[file.m$Category == "Both down",]
# head(filteredList)
# outliers1 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.1", 6], plot = FALSE)$out
# outliers2 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.2", 6], plot = FALSE)$out
# head(outliers1)
# head(outliers2)
# filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),]
# 
# c = ggplot(filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),], aes(x=factor(Category), y=value, fill = variable)) + 
#         geom_boxplot(outlier.shape = NA, width=0.05) +
#         geom_violin() +
#         geom_hline(yintercept = 0, colour="gray", linetype="dashed") +
#         ylim(-3, +3) +
#         theme_classic() 
# rm(list = c("filteredList", "outliers1", "outliers2"))
# 
# filteredList = file.m[file.m$Category == "NeoF",]
# head(filteredList)
# outliers1 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.1", 6], plot = FALSE)$out
# outliers2 = boxplot(filteredList[filteredList$variable == "median.FC.Diapause.2", 6], plot = FALSE)$out
# head(outliers1)
# head(outliers2)
# filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),]
# 
# d = ggplot(filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),], aes(x=factor(Category), y=value, fill = variable)) + 
#         geom_boxplot(outlier.shape = NA, width=0.05) +
#         geom_violin() +
#         geom_hline(yintercept = 0, colour="gray", linetype="dashed") +
#         ylim(-3, +3) +
#         theme_classic() 
# d
# rm(list = c("filteredList", "outliers1", "outliers2"))
# 
# # filteredList = file.m4[file.m3$Category == "NeoF",]
# # outliers1 = boxplot(filteredList[filteredList$variable == "median.Development.1", 6], plot = FALSE)$out
# # outliers2 = boxplot(filteredList[filteredList$variable == "median.Diapause.1", 6], plot = FALSE)$out
# # outliers3 = boxplot(filteredList[filteredList$variable == "median.Development.2", 6], plot = FALSE)$out
# # outliers4 = boxplot(filteredList[filteredList$variable == "median.Diapause.2", 6], plot = FALSE)$out
# # ggplot(filteredList[(filteredList$value %nin% outliers1 & filteredList$value %nin% outliers2),], aes(x=factor(Category), y=value, fill = variable)) + 
# #         geom_boxplot(outlier.shape = NA) + 
# #         theme_classic() 
# 
# pdf(file = "violin.pdf", width = 20)
# ggarrange(a, b, c, d + rremove("x.text"), 
#           labels = c("A", "B", "C", "D"),
#           ncol = 4, nrow = 1)
# dev.off()
# 
# 
