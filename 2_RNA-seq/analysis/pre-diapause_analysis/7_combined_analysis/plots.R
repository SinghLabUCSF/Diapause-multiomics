library("ggpubr") # To combine panels

# Set wd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library("reshape2")
library("ggplot2")

file = read.table("nfurzeri_ortho_paralogs_with_all_species.txt", header = T, sep = "\t")

pdf(file ="development_all.pdf", height = 5, width = 8)
par(mfrow=c(1,6))
# Drer paralgs ---------------------------------------------------------------------------------------------
boxplot(file$drer.10_16hpf.1[file$nfur.category == 'NeoF'],
        file$drer.10_16hpf.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "Zebrafish",
        ylab = "Expression level",
        las=1
)
text(x = 1:2,
     y = par("usr")[3] - 0.45,      ## Move labels to just below bottom of chart.
     labels = c("Gene1","Gene2"),   ## Use names from the data list.
     xpd = NA,                      ## Change the clipping region.
     srt = 45,                      ## Rotate the labels by 45 degrees.
     adj = 1.2,                     ## Adjust the labels to almost 100% right-justified.
     cex = 1.3)                     ## Increase label size.


# Olat paralgs ---------------------------------------------------------------------------------------------
boxplot(file$olat.st16.1[file$nfur.category == 'NeoF'],
        file$olat.st16.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "Medaka",
        ylab = "Expression level",
        las=1)
text(x = 1:2, y = par("usr")[3] - 0.45, labels = c("Gene1","Gene2"), xpd = NA, srt = 45, adj = 1.2, cex = 1.3)


# Alim paralgs ---------------------------------------------------------------------------------------------
boxplot(file$alim.kv30.1[file$nfur.category == 'NeoF'],
        file$alim.kv30.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "South American killifish",
        ylab = "Expression level",
        las=1)
text(x = 1:2, y = par("usr")[3] - 0.45, labels = c("Gene1","Gene2"), xpd = NA, srt = 45, adj = 1.2, cex = 1.3)


# Ast paralgs ---------------------------------------------------------------------------------------------
boxplot(file$ast.kv.1[file$nfur.category == 'NeoF'],
        file$ast.kv.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "Red-striped killifish",
        ylab = "Expression level",
        las=1)
text(x = 1:2, y = par("usr")[3] - 0.45, labels = c("Gene1","Gene2"), xpd = NA, srt = 45, adj = 1.2, cex = 1.3)


# Aaul paralgs ---------------------------------------------------------------------------------------------
boxplot(file$aaul.kv.1[file$nfur.category == 'NeoF'],
        file$aaul.kv.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "Orange killifish",
        ylab = "Expression level",
        las=1)
text(x = 1:2, y = par("usr")[3] - 0.45, labels = c("Gene1","Gene2"), xpd = NA, srt = 45, adj = 1.2, cex = 1.3)

# Nfur ---------------------------------------------------------------------------------------------
head(file)
boxplot(file$nfur.kv.1[file$nfur.category == 'NeoF'],
        file$nfur.kv.2[file$nfur.category == 'NeoF'],
        outline = F,
        main = "African turquoise killifish",
        ylab = "Expression level",
        las=1
)
text(x = 1:2, y = par("usr")[3] - 0.45, labels = c("Gene1","Gene2"), xpd = NA, srt = 45, adj = 1.2, cex = 1.3)

dev.off()

# Stats ---------------------------------------------------------------------------------------------
ks.test(file$aaul.kv.1[file$nfur.category == 'NeoF'],file$aaul.kv.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$aaul.kv.1[file$nfur.category == 'NeoF'],file$aaul.kv.2[file$nfur.category == 'NeoF'], alternative = "g")

ks.test(file$nfur.kv.1[file$nfur.category == 'NeoF'],file$nfur.kv.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$nfur.kv.1[file$nfur.category == 'NeoF'],file$nfur.kv.2[file$nfur.category == 'NeoF'], alternative = "g")

ks.test(file$ast.kv.1[file$nfur.category == 'NeoF'],file$ast.kv.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$ast.kv.1[file$nfur.category == 'NeoF'],file$ast.kv.2[file$nfur.category == 'NeoF'], alternative = "g")

ks.test(file$olat.st16.1[file$nfur.category == 'NeoF'],file$olat.st16.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$olat.st16.1[file$nfur.category == 'NeoF'],file$olat.st16.2[file$nfur.category == 'NeoF'], alternative = "g")

ks.test(file$alim.kv30.1[file$nfur.category == 'NeoF'],file$alim.kv30.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$alim.kv30.1[file$nfur.category == 'NeoF'],file$alim.kv30.2[file$nfur.category == 'NeoF'], alternative = "g")

ks.test(file$drer.10_16hpf.1[file$nfur.category == 'NeoF'],file$drer.10_16hpf.2[file$nfur.category == 'NeoF'], alternative = "l")
wilcox.test(file$drer.10_16hpf.1[file$nfur.category == 'NeoF'],file$drer.10_16hpf.2[file$nfur.category == 'NeoF'], alternative = "g")


# # All paralgs ---------------------------------------------------------------------------------------------
# boxplot(file$drer.10_16hpf.1[file$nfur.category == 'NeoF'],
#         file$drer.10_16hpf.2[file$nfur.category == 'NeoF'],
#         
#         file$olat.st16.1[file$nfur.category == 'NeoF'],
#         file$olat.st16.2[file$nfur.category == 'NeoF'],
#         
#         file$alim.kv30.1[file$nfur.category == 'NeoF'],
#         file$alim.kv30.2[file$nfur.category == 'NeoF'],
#         
#         file$ast.kv.1[file$nfur.category == 'NeoF'],
#         file$ast.kv.2[file$nfur.category == 'NeoF'],
#         
#         file$aaul.kv.1[file$nfur.category == 'NeoF'],
#         file$aaul.kv.2[file$nfur.category == 'NeoF'],
#         
#         file$nfur.kv.1[file$nfur.category == 'NeoF'],
#         file$nfur.kv.2[file$nfur.category == 'NeoF'],
#         outline = F)

