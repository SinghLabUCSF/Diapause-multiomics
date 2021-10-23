#######################################################################################
#Example Command Line Input
#Rscript Figure_3B.R
#######################################################################################

#######################################################################################
### Set enviromental variables
#######################################################################################
set.seed(1234)
CURRENT_DIR <- './'
setwd(paste(CURRENT_DIR,'Data/', sep=''))

#######################################################################################
### Load needed packages
#######################################################################################
library(ggplot2)
library(ggfortify)
library(DESeq2)

#######################################################################################
### Read in RPKM data
#######################################################################################
Input <- read.table('./RPKM_Master_0.0.txt', header = T)
Input[1] <- NULL

#######################################################################################
### Convert =data to integers and save as a separate table
#######################################################################################
for(j in 1:ncol(Input)){
  Input[,j] <- as.integer(Input[,j])
}
write.table(Input,file="INT_RPKM.txt")

#######################################################################################
### Extract species-specific comparison groups
#######################################################################################
Nfur <- cbind(Input[,1:4], Input[,6:10])
Afri <- cbind(Input[,1:4], Input[,6:10], Input[,12:16], Input[,19:20])
Kili <- cbind(Input[,1:4], Input[,6:10], Input[,12:16], Input[,19:20], Input[,21:23])
Fall <- cbind(Input[,1:4], Input[,6:10], Input[,12:16], Input[,19:20], Input[,21:23], Input[,24:27], Input[,28:31])
Dia <- cbind(Input[,1:4], Input[,6:10],Input[,21:23])

#######################################################################################
### Remove values missing within at least one species
#######################################################################################
Nfur <- na.omit(Nfur)
Afri <- na.omit(Afri)
Kili <- na.omit(Kili)
Fall <- na.omit(Fall)
Dia <- na.omit(Dia)

#######################################################################################
### Convert all data types to 'numeric'
#######################################################################################
for (i in 1:ncol(Afri)){
  Afri[,i] <- as.numeric(Afri[,i])
}
for (i in 1:ncol(Kili)){
  Kili[,i] <- as.numeric(Kili[,i])
}
for (i in 1:ncol(Fall)){
  Fall[,i] <- as.numeric(Fall[,i])
}
for (i in 1:ncol(Dia)){
  Dia[,i] <- as.numeric(Dia[,i])
}

#######################################################################################
### Apply VST normalization to data
#######################################################################################
Nfur <- vst(as.matrix(Nfur))
Afri <- vst(as.matrix(Afri))
Kili <- vst(as.matrix(Kili))
Fall <- varianceStabilizingTransformation(as.matrix(Fall))
Dia <- vst(as.matrix(Dia))

#######################################################################################
### Transpose table for use by prcomp
#######################################################################################
Nfur <- t(Nfur)
Afri <- t(Afri)
Kili <- t(Kili)
Fall <- t(Fall)
Dia <- t(Dia)

counter <- c(ncol(Nfur), ncol(Afri), ncol(Kili), ncol(Fall), ncol(Dia))

#######################################################################################
### Generate eigan vectors for PCA
#######################################################################################
test1 <- prcomp(Nfur, scale.=T)
test2 <- prcomp(Afri, scale.=T)
test3 <- prcomp(Kili, scale.=T)
test4 <- prcomp(Fall, scale.=T)
test5 <- prcomp(Dia, scale.=T)

#######################################################################################
### Label all samples for graphing
#######################################################################################
nNfur <- c(rep('Nfur', 9))
nAaus <- c(rep('Aaus', 4))
nAstr <- c(rep('Astr', 3))
nAlim <- c(rep('Alim', 3))
nOlat <- c(rep('Olat', 4))
nDrer <- c(rep('Drer', 4))
nAfri <- c(nNfur, nAaus, nAstr)
nKili <- c(nNfur, nAaus, nAstr, nAlim)
nFall <- c(nNfur, nAaus, nAstr, nAlim, nOlat, nDrer)
nDia <- c(nNfur,nAlim)

sNfur <- c('Nfur_Dia6', 'Nf_Dia6', 'Nf_esc', 'Nf_esc', 'Nf_KvY', 'Nf_KvY', 'Nf_DiaM', 'Nf_DiaM', 'Nf_DiaM')
sAaus <- c('Aa_Kv', 'Aa_Kv', 'Aa_HB', 'Aa_Kv')
sAstr <- c('As_HB', 'As_KvY', 'As_KvY')
sAlim <- c('Al_DiaM', 'Al_esc', 'Al_Kv')
sOlat <- c('Ol_19', 'Ol_19', 'Ol_25', 'Ol_25')
sDrer <- c('Dr_8s', 'Dr_8s','Dr_48hr', 'Dr_48hr')
sAfri <- c(sNfur, sAaus, sAstr)
sKili <- c(sNfur, sAaus, sAstr, sAlim)
sFall <- c(sNfur, sAaus, sAstr, sAlim, sOlat, sDrer)
sDia <- c(sNfur,sAlim)

dNfur <- c('Dia', 'Dia', 'Dev', 'Dev', 'Dev', 'Dev', 'Dia', 'Dia', 'Dia')
dAaus <- c('Dev', 'Dev', 'Dev', 'Dev')
dAstr <- c('Dev', 'Dev', 'Dev')
dAlim <- c('Dia', 'Dev', 'Dev')
dOlat <- c('Dev', 'Dev', 'Dev', 'Dev')
dDrer <- c('Dev', 'Dev', 'Dev', 'Dev')
dAfri <- c(dNfur, dAaus, dAstr)
dKili <- c(dNfur, dAaus, dAstr, dAlim)
dFall <- c(dNfur, dAaus, dAstr, dAlim, dOlat, dDrer)
dDia <- c(dNfur,dAlim)

Nfur <- cbind(nNfur, sNfur, dNfur, Nfur)
Afri <- cbind(nAfri, sAfri, dAfri, Afri)
Kili <- cbind(nKili, sKili, dKili, Kili)
Fall <- cbind(nFall, sFall, dFall, Fall)
Dia <- cbind(nDia, sDia, dDia, Dia)

#######################################################################################
### Extract Diapause/Development PC2 loadings and save as a table
#######################################################################################
last <- as.data.frame(test5$rotation[,2])
write.table(last, './Nfur_Alim_PC2.txt', sep='\t')

#######################################################################################
### Graph and save PCAs
#######################################################################################
pdf('Figure3_B1.pdf')
autoplot(test1, x=1, y=2, data=Nfur, size=10, colour='dNfur', shape='nNfur') +
  scale_colour_manual(values=c('orange', 'light blue')) +
  scale_shape_manual(values=c(19)) +
  labs(title= 'N. furzeri PCA')
dev.off()

pdf('Figure3_B2.pdf')
autoplot(test3, x=1, y=2, data=Kili, size=10, colour='dKili', shape='nKili') +
  scale_colour_manual(values=c('orange', 'light blue')) +
  scale_shape_manual(values=c(18, 15, 17, 19)) +
  labs(title= 'Killifish PCA')
dev.off()

pdf('Figure3_B3.pdf')
autoplot(test5, x=1, y=2, data=Dia, size=10, colour='dDia', shape='nDia') +
  scale_colour_manual(values=c('orange', 'light blue')) +
  scale_shape_manual(values=c(15, 19)) +
  labs(title= 'Diapause-Capable Killifish PCA')
dev.off()










