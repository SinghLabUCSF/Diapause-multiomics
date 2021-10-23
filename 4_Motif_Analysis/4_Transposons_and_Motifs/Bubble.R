#######################################################################################
#Example Command Line Input
#Rscript Bubble.R
#######################################################################################

#######################################################################################
### Set enviromental variables
#######################################################################################
set.seed(1234)

CURRENT_DIR <- './'

setwd(CURRENT_DIR)

library('stats')
library('MutationalPatterns')
library('ggplot2')
library('circlize')

#######################################################################################
### Read in data and adjust headers with special characters
#######################################################################################
Input <- read.table(paste(CURRENT_DIR, 'Data/TE_Family_Data_Filtered.txt', sep=''), header = T)
Input[29,1] <- 'Piggybac-like element'
Input[42,1] <- 'hAT-like element'
Input[44,1] <- 'CMC-EnSpm-like element'
Input[48,1] <- 'MULE-MuDR-like element'
Input[58,1] <- 'PIF-Harbinger-like element'
Input[60,1] <- 'Zisupton-like element'
Input[66,1] <- 'hAT-Ac-like element'
Input[67,1] <- 'TcMar-Stowaway-like element'
Input[68,1] <- 'Kolobok-T2-like element'
Input[72,1] <- 'SINE-like element'
Input[77,1] <- 'hAT-Tip100-like element'
Input[78,1] <- 'Zator-like element'
Input[81,1] <- 'hAT-Charlie-like element'
Input[83,1] <- 'Kolobok-like element'
Input[84,1] <- 'TcMar-like element'
Input[72,10] <- 'SINE'

#######################################################################################
### Generate empty structures for data and populate each of them with data subsets
#######################################################################################
TOT <- Input[nrow(Input),6]
Summary <- NULL
Namer <- NULL
TEer <- NULL
DEer <- NULL
Typer <- NULL
Clas <- NULL

for (i in 1:nrow(Input)){
  temp1 <- binomial_test(Input[i,3], TOT, Input[i,6])
  if (temp1[1] == 'enrichment'){
    temp1 <- cbind(temp1[2],(Input[i,7]/Input[i,3]))
    colnames(temp1) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary, temp1)
  } else if (temp1[1] != 'enrichment'){
    temp1 <- cbind(temp1[2],-1*(Input[i,3]/Input[i,7]))
    #temp1 <- cbind(1,0)
    colnames(temp1) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary, temp1)
  }
  
  temp2 <- binomial_test(Input[i,5], TOT, Input[i,6])
  if (temp2[1] == 'enrichment'){
    temp2 <- cbind(temp2[2],(Input[i,7]/Input[i,5]))
    colnames(temp2) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary,temp2)
  } else if (temp2[1] != 'enrichment'){
    temp2 <- cbind(temp2[2],-1*(Input[i,5]/Input[i,7]))
    #temp2 <- cbind(1,0)
    colnames(temp2) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary,temp2)
  }
  
  temp3 <- binomial_test(Input[i,9], TOT, Input[i,6])
  if (temp3[1] == 'enrichment'){
    temp3 <- cbind(temp3[2],(Input[i,7]/Input[i,9]))
    colnames(temp3) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary,temp3)
  } else if (temp3[1] != 'enrichment'){
    temp3 <- cbind(temp3[2],-1*(Input[i,9]/Input[i,7]))
    #temp3 <- cbind(1,0)
    colnames(temp3) <- c('Pval', 'FoldC')
    Summary <- rbind(Summary,temp3)
  }
} 

#######################################################################################
### Summarize data and stats while correcting P-values
#######################################################################################
for (i in 1:nrow(Summary)){
  if (is.na(Summary[i,2])){
    Summary[i,2] <- 0
  } else if (is.infinite(Summary[i,2]) & Summary[i,2]*-1 < 0){
    Summary[i,2] <- 7
  } else if (is.infinite(Summary[i,2]) & Summary[i,2]*-1 > 0){
    Summary[i,2] <- -7
  }
}

temp1 <- p.adjust(Summary[,1], method='BH')
Summary<- cbind(Summary[,1], temp1, Summary[,2])


#######################################################################################
### Recombine all data for graphing
#######################################################################################
for (i in 1:nrow(Input)){
    for (j in 1:3){
      Namer <- rbind(Namer,i)
      TEer <- rbind(TEer,Input[i,1])
      DEer <- rbind(DEer,Input[i,6])
      Clas <- rbind(Clas,Input[i,10])
    }
    Typer <- rbind(Typer, 'Genomic abundance')
    Typer <- rbind(Typer, 'Chromatin abundance')
    Typer <- rbind(Typer, 'Control loci abundance')
}

Summary <- cbind(TEer,Clas,Namer,Typer,DEer,Summary)
colnames(Summary) <- c('TE','Class','TEN','Test','Count','Pval','Padj','FC')
Summary <- as.data.frame(Summary)


#######################################################################################
### Keep TE if any condition is significant after multiple test correction
#######################################################################################
clean <- NULL
counter <- 0
for(i in 1:nrow(Summary)){
  if (Summary[i,7] < 0.1 & counter == 0){
    if (i%%3 == 1){
      clean <- rbind(clean, Summary[i,],Summary[i+1,], Summary[i+2,])
      counter <- 2
    } else if(i%%3 == 2){
      clean <- rbind(clean, Summary[i-1,], Summary[i,], Summary[i+1,])
      counter <- 1
    } else if (i%%3 == 0){
      clean <- rbind(clean, Summary[i-2,], Summary[i-1,], Summary[i,])
      counter <- 0
    }
  } else if (counter != 0){
    counter <- counter - 1
  }
}

Namer2 <- NULL
for (i in 1:nrow(clean)){
  Namer2 <- rbind(Namer2,i)
}

clean <- cbind(clean[,1:2],Namer2,clean[,4:8])

for (i in 1:nrow(clean)){
  if (clean[i,7] == 0){
    clean[i,7] <- 1.0e-35
  }
}


#######################################################################################
### Subset and Label data for graphing
#######################################################################################
clean$Padj <- as.numeric(clean$Padj)
clean$FC <- as.numeric(clean$FC)


T1 <- subset(clean, subset=(Test == 'Genomic abundance'))
T2 <- subset(clean, subset=(Test == 'Chromatin abundance'))
T3 <- subset(clean, subset=(Test == 'Control loci abundance'))

T1S <- subset(clean, subset=(Test == 'Genomic abundance' & as.numeric(Padj) < 0.1))
T2S <- subset(clean, subset=(Test == 'Chromatin abundance' & as.numeric(Padj) < 0.1))
T3S <- subset(clean, subset=(Test == 'Control loci abundance' & as.numeric(Padj) < 0.1))


#######################################################################################
### Color Data by fold change and P-value
#######################################################################################
T2C <- NULL
for (i in 1:nrow(T2S)){
  if (T2S[i,1] %in% T1S$TE & T2S[i,1] %in% T3S$TE){
    T2C <- rbind(T2C, T2S[i,])
  }
}

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T1)){
  Namer <- rbind(Namer, i)
  if (T1[i,8] < 0 & T1[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T1[i,8] > 0 & T1[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T1$Namer <- Namer
T1$Paint <- Paint

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T2)){
  Namer <- rbind(Namer, i)
  if (T2[i,8] < 0 & T2[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T2[i,8] > 0 & T2[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T2$Namer <- Namer
T2$Paint <- Paint

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T3)){
  Namer <- rbind(Namer, i)
  if (T3[i,8] < 0 & T3[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T3[i,8] > 0 & T3[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T3$Namer <- Namer
T3$Paint <- Paint

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T1S)){
  Namer <- rbind(Namer, i)
  if (T1S[i,8] < 0 & T1S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T1S[i,8] > 0 & T1S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T1S$Namer <- Namer
T1S$Paint <- Paint

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T2S)){
  Namer <- rbind(Namer, i)
  if (T2S[i,8] < 0 & T2S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T2S[i,8] > 0 & T2S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T2S$Namer <- Namer
T2S$Paint <- Paint

Namer <-NULL
Paint <- NULL
for (i in 1:nrow(T3S)){
  Namer <- rbind(Namer, i)
  if (T3S[i,8] < 0 & T3S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Blue')
  }else if (T3S[i,8] > 0 & T3S[i,7] <= 0.1){
    Paint <- rbind(Paint, 'Red')
  } else {
    Paint <- rbind(Paint, 'Grey')
  }
}
T3S$Namer <- Namer
T3S$Paint <- Paint

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T1dna <- subset(T1, subset=(Class == 'DNA'))
T1sine <- subset(T1, subset=(Class == 'SINE'))
T1line <- subset(T1, subset=(Class == 'LINE'))

T1dna <- T1dna[order(T1dna$FC),]
T1sine <- T1sine[order(T1sine$FC),]
T1line <- T1line[order(T1line$FC),]

Newsine <- seq(1,nrow(T1sine),1)
Newline <- seq((1+nrow(T1sine)),(nrow(T1sine)+nrow(T1line)),1)
Newdna <- seq((1+nrow(T1sine)+nrow(T1line)),(nrow(T1sine)+nrow(T1line)+nrow(T1dna)),1)
Neword <- c(Newsine,Newline,Newdna)

T1final <- rbind(T1sine,T1line,T1dna)
T1final$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_1All.pdf', sep=''))
ggplot(data=T1final, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T1final), 1)) +
  scale_size(range=c(1,6)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=1.5) +
  geom_hline(yintercept=9.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=36), fill=NA, color='Black', size=0.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=36), fill='White', color=NA, size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue", "grey", "red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T1Sdna <- subset(T1S, subset=(Class == 'DNA'))
T1Ssine <- subset(T1S, subset=(Class == 'SINE'))
T1Sline <- subset(T1S, subset=(Class == 'LINE'))

T1Sdna <- T1Sdna[order(T1Sdna$FC),]
T1Ssine <- T1Ssine[order(T1Ssine$FC),]
T1Sline <- T1Sline[order(T1Sline$FC),]

Newline <- seq(1,nrow(T1Sline),1)
Newdna <- seq((1+nrow(T1Sline)),(nrow(T1Sline)+nrow(T1Sdna)),1)
Neword <- c(Newline,Newdna)

T1Sfinal <- rbind(T1Sline,T1Sdna)
T1Sfinal$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_1Sig.pdf', sep=''))
ggplot(data=T1Sfinal, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T1Sfinal), 1)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0)) +
  scale_size(range=c(1,10)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=6.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=17), fill='White', color=NA, size=0.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=17), fill=NA, color='Black', size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue","red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T2dna <- subset(T2, subset=(Class == 'DNA'))
T2sine <- subset(T2, subset=(Class == 'SINE'))
T2line <- subset(T2, subset=(Class == 'LINE'))

T2dna <- T2dna[order(T2dna$FC),]
T2sine <- T2sine[order(T2sine$FC),]
T2line <- T2line[order(T2line$FC),]

Newsine <- seq(1,nrow(T2sine),1)
Newline <- seq((1+nrow(T2sine)),(nrow(T2sine)+nrow(T2line)),1)
Newdna <- seq((1+nrow(T2sine)+nrow(T2line)),(nrow(T2sine)+nrow(T2line)+nrow(T2dna)),1)
Neword <- c(Newsine,Newline,Newdna)

T2final <- rbind(T2sine,T2line,T2dna)
T2final$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_2All.pdf', sep=''))
ggplot(data=T2final, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T2final), 1)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0))+
  scale_size(range=c(1,6)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=1.5) +
  geom_hline(yintercept=9.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=36), fill='White', color=NA, size=0.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=36), fill=NA, color='Black', size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue", "grey", "red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T2Sdna <- subset(T2S, subset=(Class == 'DNA'))
T2Ssine <- subset(T2S, subset=(Class == 'SINE'))
T2Sline <- subset(T2S, subset=(Class == 'LINE'))

T2Sdna <- T2Sdna[order(T2Sdna$FC),]
T2Ssine <- T2Ssine[order(T2Ssine$FC),]
T2Sline <- T2Sline[order(T2Sline$FC),]

Newline <- seq(1,nrow(T2Sline),1)
Newdna <- seq((1+nrow(T2Sline)),(nrow(T2Sline)+nrow(T2Sdna)),1)
Neword <- c(Newline,Newdna)

T2Sfinal <- rbind(T2Sline,T2Sdna)
T2Sfinal$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_2Sig.pdf', sep=''))
ggplot(data=T2Sfinal, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T2Sfinal), 1)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0))+
  scale_size(range=c(1,10)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=3.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=12), fill=NA, color='Black', size=0.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=12), fill='White', color=NA, size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue","red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T3dna <- subset(T3, subset=(Class == 'DNA'))
T3sine <- subset(T3, subset=(Class == 'SINE'))
T3line <- subset(T3, subset=(Class == 'LINE'))

T3dna <- T3dna[order(T3dna$FC),]
T3sine <- T3sine[order(T3sine$FC),]
T3line <- T3line[order(T3line$FC),]

Newsine <- seq(1,nrow(T3sine),1)
Newline <- seq((1+nrow(T3sine)),(nrow(T3sine)+nrow(T3line)),1)
Newdna <- seq((1+nrow(T3sine)+nrow(T3line)),(nrow(T3sine)+nrow(T3line)+nrow(T3dna)),1)
Neword <- c(Newsine,Newline,Newdna)

T3final <- rbind(T3sine,T3line,T3dna)
T3final$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_3All.pdf', sep=''))
ggplot(data=T3final, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T3final), 1)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0))+
  scale_size(range=c(1,6)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=1.5) +
  geom_hline(yintercept=9.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=36), fill='White', color=NA, size=0.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=36), fill=NA, color='Black', size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue", "grey", "red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
T3Sdna <- subset(T3S, subset=(Class == 'DNA'))
T3Ssine <- subset(T3S, subset=(Class == 'SINE'))
T3Sline <- subset(T3S, subset=(Class == 'LINE'))

T3Sdna <- T3Sdna[order(T3Sdna$FC),]
T3Ssine <- T3Ssine[order(T3Ssine$FC),]
T3Sline <- T3Sline[order(T3Sline$FC),]

Newsine <- seq(1,nrow(T3Ssine),1)
Newline <- seq((1+nrow(T3Ssine)),(nrow(T3Ssine)+nrow(T3Sline)),1)
Newdna <- seq((1+nrow(T3Ssine)+nrow(T3Sline)),(nrow(T3Ssine)+nrow(T3Sline)+nrow(T3Sdna)),1)
Neword <- c(Newsine,Newline,Newdna)

T3Sfinal <- rbind(T3Ssine,T3Sline,T3Sdna)
T3Sfinal$Namer2 <- Neword

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_3Sig.pdf', sep=''))
ggplot(data=T3Sfinal, aes(x=as.numeric(FC), y=as.numeric(Namer2), size=as.numeric(-log(Padj)))) +
  geom_point(aes(color=Paint)) +
  scale_y_continuous(expand=c(0,0), minor_breaks = seq(0, nrow(T3final), 1)) +
  scale_x_continuous(limits=c(-16,8), breaks=c(-8,-6,-4,-2,0,2,4,6,8), expand=c(0,0)) +
  scale_size(range=c(1,10)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=1.5) +
  geom_hline(yintercept=8.5) +
  geom_rect(mapping=aes(xmin=-16, xmax=-8, ymin=0, ymax=21), fill='White', color=NA, size=0.5) +
  geom_rect(mapping=aes(xmin=-8, xmax=8, ymin=0, ymax=21), fill=NA, color='Black', size=0.5) +
  labs(title = 'Transpsoable Element Enrichment', x='Fold Change') +
  geom_text(aes(label = paste(TE, '-', sep=' '), x = -8, y = Namer2), size=3, hjust=1) +
  theme(panel.grid.major = element_line(colour="light grey", size=0.1), panel.grid.minor = element_line(size=0.1, colour="light grey"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_color_manual(values=c("blue", "red"))
dev.off()

#######################################################################################
### Reorder data by TE Super-Class
#######################################################################################
FCs <- NULL
Super <- rbind(T1final,T2final,T3final)
Super <- Super[order(Super$Class),]
Super$TE <- factor(Super$TE, levels=T1final$TE)
Super$Test <- factor(Super$Test, levels=c('Control loci abundance', 'Chromatin abundance', 'Genomic abundance'))
for (i in 1:nrow(Super)){
  if (as.numeric(Super[i,7])>0.1){
    FCs <- rbind(FCs, NA)
  } else {
    FCs <- rbind(FCs, Super[i,8])
  }
}
Super <- cbind(Super, FCs)

PAL <- colorRamp2(c(-5,-3,-1,0,1,3,5), c( '#0000FF','#6666FF','#9999CC','#CCCCCC','#CC9999','#FF6666','#FF0000'))
TPAL <- PAL(c(-5,-3,-1,0,1,3,5,7,9,11,13,15))
pdf(paste(CURRENT_DIR, 'Data/TE_Bub_Pile.pdf', sep = ''), height = 12,width = 5)
ggplot(Super, aes(x = Test, y = TE, color=FC)) +
  geom_point(aes(size = -log10(Padj))) + 
  scale_size_continuous(breaks=c(0,1,10,20,30), labels=c('P=1','P=0.1','P=1e-10','P=1e-20','P=1e-30')) +
  scale_color_gradientn(colors=c(TPAL)) +  
  theme_bw() +
  theme(axis.title= element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
dev.off()


Summary$Padj <- as.numeric(Summary$Padj)
Summary$FC <- as.numeric(Summary$FC)
piled <- subset(Summary, subset=(Class == 'DNA'))
piled <- subset(piled, subset=(Test == 'Genomic abundance'))
piles <- subset(Summary, subset=(Class == 'SINE'))
piles <- subset(piles, subset=(Test == 'Genomic abundance'))
pilel <- subset(Summary, subset=(Class == 'LINE'))
pilel <- subset(pilel, subset=(Test == 'Genomic abundance'))

pileg <- rbind(piles,pilel,piled)
Summary$TE <- factor(Summary$TE, levels=pileg$TE)
for (i in 1:nrow(Summary)){
  if (Summary[i,7] == 0){
    Summary[i,7] <- 1.0e-35
  }
}

SSummary <- NULL
for (i in 1:nrow(Summary)){
  if (Summary[i,2] != 'Class'){
    SSummary <- rbind(SSummary, Summary[i,])
  }
}


FCs <- NULL
for (i in 1:nrow(SSummary)){
  if (as.numeric(SSummary[i,7])>0.1){
    FCs <- rbind(FCs, NA)
  } else {
    FCs <- rbind(FCs, SSummary[i,8])
  }
}
SSummary <- cbind(SSummary, FCs)
SSummary$Test <- factor(SSummary$Test, levels=c('Control loci abundance', 'Chromatin abundance', 'Genomic abundance'))

#######################################################################################
### Graph Data
#######################################################################################
pdf(paste(CURRENT_DIR,'Data/TE_Bub_Pile_All.pdf', spe=''), height = 18 ,width = 5)
ggplot(SSummary, aes(x =Test, y = TE, color=FCs)) +
  scale_color_gradientn(colors=c(TPAL)) + 
  scale_size_continuous(breaks=c(0,1,10,20,30), labels=c('P=1','P=0.1','P=1e-10','P=1e-20','P=1e-30')) +
  geom_point(aes(size = -log10(Padj))) + 
  theme_bw() +
  theme(axis.title= element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
dev.off()

