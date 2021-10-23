#######################################################################################
#Example Command Line Input
#Rscript Figure_3C_S10.R
#######################################################################################

#######################################################################################
### Set enviromental variables
#######################################################################################
set.seed(1234)
CURRENT_DIR <- './'
USER <- paste(CURRENT_DIR,'Data/', sep = '')
setwd(USER)

#######################################################################################
### Load needed packages
#######################################################################################
library('ggplot2')
library('dplyr')
library('tidyverse')

#######################################################################################
### Define function
#######################################################################################
# This function counts the frequency of give groups withing the data and calculates
# that number as a precentage
#######################################################################################

Runner <- function(IN_D,NAMER) {
  ct = count(IN_D, IN_D$Con)
  colnames(ct) = c("Con", "freq")
  ct$percent = round((ct$freq*100/sum(ct$freq)), digits = 1)
  ct$joint = NAMER
  total = 100
  for(i in 1:nrow(ct)){
    j = total - ct[i,3]/2
    ct$position[i] = j
    print (j)
    total = total - ct[i,3]
  }
  return(ct)
}

#######################################################################################
### Loop through both alignment and peak data
#######################################################################################
TYPER <- c('AS', 'PS')

for (i in 1:length(TYPER)){
  SHOUT <- TYPER[i]
  
  #######################################################################################
  ### Define loop-specific variables
  #######################################################################################
  if(SHOUT == 'AS'){
    SETTER <- 'Alignment'
    orders1 <- c('Exon', 'UTR', 'Promoter', 'Intron', 'Intergenic')
    OUTER <- 'Alignment_Simplified'
    Plate <- c("#CC0A0A", "#FC7C0C","#0ACC81")
    Door <- './Final_Con_List_align_Edit_R.txt'
    
  }else if(SHOUT == 'PS'){
    SETTER <- 'Peak'
    orders1 <- c('Promoter', 'Exon', 'UTR', 'Intron', 'Intergenic')
    OUTER <- 'Peak_Simplified'
    Plate <- c("#CC0A0A", "#FC7C0C","#0ACC81")
    Door <- './Final_Con_List_peak_Edit_R.txt'
  }
  
  #######################################################################################
  ### Read in data and separate it by various metadata statuses
  #######################################################################################
  Input <- read.table(Door, header = T)
  
  DATA <- as.data.frame(Input)
  DE <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia'])
  NDE <- data.frame(Con=DATA$Con[DATA$DE != 'DE_Dia'])
  NeoF <- data.frame(Con=DATA$Con[DATA$Stat == 'N1' | DATA$Stat =='N2'])
  NNeoF <- data.frame(Con=DATA$Con[DATA$Stat != 'N1' & DATA$Stat !='N2'])
  ANC <- data.frame(Con=DATA$Con[DATA$age == 'V'])
  OLD <- data.frame(Con=DATA$Con[DATA$age == 'F'])
  REC <- data.frame(Con=DATA$Con[DATA$age == 'K'])
  PRO <- data.frame(Con=DATA$Con[DATA$Type == 'Promoter'])
  EXN <- data.frame(Con=DATA$Con[DATA$Type == 'Exon'])
  INN <- data.frame(Con=DATA$Con[DATA$Type == 'Intron'])
  UTR <- data.frame(Con=DATA$Con[DATA$Type == 'UTR'])
  IGC <- data.frame(Con=DATA$Con[DATA$Type == 'Intergenic'])
  
  #######################################################################################
  ### Further sub-divide data by metadata statuses
  #######################################################################################
  
  NeoF_DE <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$DE == 'DE_Dia')])
  NNeoF_DE <- data.frame(Con=DATA$Con[(DATA$Stat != 'N1' & DATA$Stat !='N2') & (DATA$DE == 'DE_Dia')])
  
  NeoF_PRO <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Promoter')])
  NeoF_EXN <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Exon')])
  NeoF_INN <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Intron')])
  NeoF_UTR <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'UTR')])
  NeoF_IGC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Intergenic')])
  
  NeoF_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'V')])
  NeoF_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'F')])
  NeoF_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'K')])
  
  DE_PRO  <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$Type == 'Promoter'])
  DE_EXN  <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$Type == 'Exon'])
  DE_INN  <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$Type == 'Intron'])
  DE_UTR  <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$Type == 'UTR'])
  DE_IGC  <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$Type == 'Intergenic'])
  
  DE_ANC <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$age == 'V'])
  DE_OLD <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$age == 'F'])
  DE_REC <- data.frame(Con=DATA$Con[DATA$DE == 'DE_Dia' & DATA$age == 'K'])
  
  #######################################################################################
  
  NeoF_DE_PRO <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Promoter') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_EXN <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Exon') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_INN <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Intron') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_UTR <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'UTR') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_IGC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$Type == 'Intergenic') & (DATA$DE == 'DE_Dia')])
  
  NeoF_DE_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'V') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'F') & (DATA$DE == 'DE_Dia')])
  NeoF_DE_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1' | DATA$Stat =='N2') & (DATA$age == 'K') & (DATA$DE == 'DE_Dia')])
  
  #######################################################################################
  
  Neo1_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'V')])
  Neo1_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'F')])
  Neo1_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'K')])
  
  Neo2_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'V')])
  Neo2_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'F')])
  Neo2_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'K')])
  
  Neo1_DE_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'V') & (DATA$DE == 'DE_Dia')])
  Neo1_DE_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'F') & (DATA$DE == 'DE_Dia')])
  Neo1_DE_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N1') & (DATA$age == 'K') & (DATA$DE == 'DE_Dia')])
  
  Neo2_DE_ANC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'V') & (DATA$DE == 'DE_Dia')])
  Neo2_DE_OLD <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'F') & (DATA$DE == 'DE_Dia')])
  Neo2_DE_REC <- data.frame(Con=DATA$Con[(DATA$Stat == 'N2') & (DATA$age == 'K') & (DATA$DE == 'DE_Dia')])
  
  #######################################################################################
  ### Generate precentage counts for each sub-divided data category
  #######################################################################################
  DATA2 <- Runner(DATA,'Genome-Wide')
  DE2 <- Runner(DE, 'DE-Diapause')
  NDE2 <- Runner(NDE, 'Other')
  NeoF2 <- Runner(NeoF, 'NeoF')
  NNeoF2 <- Runner(NNeoF, 'Other')
  ANC2 <- Runner(ANC, 'Vertabrate-Origin')
  OLD2 <- Runner(OLD, 'Fish-Origin')
  REC2 <- Runner(REC, 'Killifish-Origin')
  PRO2 <- Runner(PRO, 'Promoter')
  EXN2 <- Runner(EXN, 'Exon')
  INN2 <- Runner(INN, 'Intron')
  UTR2 <- Runner(UTR, 'UTR')
  IGC2 <- Runner(IGC, 'Intergenic')
  
  #######################################################################################
  
  NeoF_DE2 <- Runner(NeoF_DE, 'NeoF')
  NNeoF_DE2 <- Runner(NNeoF_DE, 'Other')
  
  NeoF_PRO2 <- Runner(NeoF_PRO, 'Promoter')
  NeoF_EXN2 <- Runner(NeoF_EXN, 'Exon')
  NeoF_INN2 <- Runner(NeoF_INN, 'Intron')
  NeoF_UTR2 <- Runner(NeoF_UTR, 'UTR')
  NeoF_IGC2 <- Runner(NeoF_IGC, 'Intergenic')
  
  NeoF_ANC2 <- Runner(NeoF_ANC, 'Vertabrate-Origin')
  NeoF_OLD2 <- Runner(NeoF_OLD, 'Fish-Origin')
  NeoF_REC2 <- Runner(NeoF_REC, 'Killifish-Origin')
  
  DE_PRO2 <- Runner(DE_PRO, 'Promoter')
  DE_EXN2 <- Runner(DE_EXN, 'Exon')
  DE_INN2 <- Runner(DE_INN, 'Intron')
  DE_UTR2 <- Runner(DE_UTR, 'UTR')
  DE_IGC2 <- Runner(DE_IGC, 'Intergenic')
  
  DE_ANC2 <- Runner(DE_ANC, 'Vertabrate-Origin')
  DE_OLD2 <- Runner(DE_OLD, 'Fish-Origin')
  DE_REC2 <- Runner(DE_REC, 'Killifish-Origin')
  
  #######################################################################################
  
  NeoF_DE_PRO2 <- Runner(NeoF_DE_PRO, 'Promoter')
  NeoF_DE_EXN2 <- Runner(NeoF_DE_EXN, 'Exon')
  NeoF_DE_INN2 <- Runner(NeoF_DE_INN, 'Intron')
  NeoF_DE_UTR2 <- Runner(NeoF_DE_UTR, 'UTR')
  NeoF_DE_IGC2 <- Runner(NeoF_DE_IGC, 'Intergenic')
  
  NeoF_DE_ANC2 <- Runner(NeoF_DE_ANC, 'Vertabrate-Origin')
  NeoF_DE_OLD2 <- Runner(NeoF_DE_OLD, 'Fish-Origin')
  NeoF_DE_REC2 <- Runner(NeoF_DE_REC, 'Killifish-Origin')
  
  #######################################################################################
  
  Neo1_ANC2 <- Runner(Neo1_ANC, 'Vertabrate-Origin (N1)')
  Neo1_OLD2 <- Runner(Neo1_OLD, 'Fish-Origin (N1)')
  Neo1_REC2 <- Runner(Neo1_REC, 'Killifish-Origin (N1)')
  
  Neo2_ANC2 <- Runner(Neo2_ANC, 'Vertabrate-Origin (N2)')
  Neo2_OLD2 <- Runner(Neo2_OLD, 'Fish-Origin (N2)')
  Neo2_REC2 <- Runner(Neo2_REC, 'Killifish-Origin (N2)')
  
  Neo1_DE_ANC2 <- Runner(Neo1_DE_ANC, 'Vertabrate-Origin (N1)')
  Neo1_DE_OLD2 <- Runner(Neo1_DE_OLD, 'Fish-Origin (N1)')
  Neo1_DE_REC2 <- Runner(Neo1_DE_REC, 'Killifish-Origin (N1)')
  
  Neo2_DE_ANC2 <- Runner(Neo2_DE_ANC, 'Vertabrate-Origin (N2)')
  Neo2_DE_OLD2 <- Runner(Neo2_DE_OLD, 'Fish-Origin (N2)')
  Neo2_DE_REC2 <- Runner(Neo2_DE_REC, 'Killifish-Origin (N2)')
  
  #######################################################################################
  ### Combine sub-divided data for graphing purposes and output as a table
  #######################################################################################
  G1 <- DATA2
  G2 <- rbind(DE2,NDE2)
  G3 <- rbind(NeoF2,NNeoF2)
  G4 <- rbind(ANC2,OLD2,REC2)
  G5 <- rbind(PRO2,EXN2,INN2,UTR2,IGC2)
  G6 <- rbind(NeoF_DE2,NNeoF_DE2)
  G7 <- rbind(NeoF_PRO2,NeoF_EXN2,NeoF_INN2,NeoF_UTR2,NeoF_IGC2)
  G8 <- rbind(NeoF_ANC2,NeoF_OLD2,NeoF_REC2)
  G9 <- rbind(DE_PRO2,DE_EXN2,DE_INN2,DE_UTR2,DE_IGC2)
  G10 <- rbind(DE_ANC2,DE_OLD2,DE_REC2)
  G11 <- rbind(NeoF_DE_PRO2,NeoF_DE_EXN2,NeoF_DE_INN2,NeoF_DE_UTR2,NeoF_DE_IGC2)
  G12 <- rbind(NeoF_DE_ANC2,NeoF_DE_OLD2,NeoF_DE_REC2)
  G13 <- rbind(DATA2,DE2,NeoF_DE2,NeoF_DE_ANC2)
  G14 <- rbind(Neo1_ANC2,Neo2_ANC2,Neo1_OLD2,Neo2_OLD2,Neo1_REC2,Neo2_REC2)
  G15 <- rbind(Neo1_DE_ANC2,Neo2_DE_ANC2,Neo1_DE_OLD2,Neo2_DE_OLD2,Neo1_DE_REC2,Neo2_DE_REC2)
  GA <- rbind(G1,G2,G3,G4,G5,G6,G9,G10,G7,G8,G11,G12,G13,G14,G15)
  GL <- c(rep('Graph1',nrow(G1)), rep('Graph2',nrow(G2)), rep('Graph3',nrow(G3)), rep('Graph4',nrow(G4)),
          rep('Graph5',nrow(G5)), rep('Graph6',nrow(G6)), rep('Graph7',nrow(G9)), rep('Graph8',nrow(G10)),
          rep('Graph9',nrow(G7)), rep('Graph10',nrow(G8)), rep('Graph11',nrow(G11)), rep('Graph12',nrow(G12)),
          rep('Graph13',nrow(G13)), rep('Graph14',nrow(G14)), rep('Graph15',nrow(G15)))
  GAE <- cbind(GA,GL)
  
  write.table(GAE, file=paste(OUTER,'_Data.txt', sep=''),sep='\t')
  
  #######################################################################################
  ### Graph data
  #######################################################################################
  
  pdf(file = paste(OUTER,'_Horizontal.pdf',sep=''))
  
  p <- ggplot(G2, aes(x=joint, y=percent, fill=Con)) +
    geom_bar(stat='identity') +
    labs(x='Accessibility Status',y='Percentage',title = paste(SETTER,' Conservation Breakdown By Accessibility', sep='')) +
    scale_fill_manual(values = Plate) +
    coord_flip() +
    scale_x_discrete(limits=c('DE-Diapause','Other')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom") +
    geom_text(aes(label = sprintf("%1.2f%%", percent), y = position))
  print(p)
  
  p <- ggplot(G6, aes(x=joint, y=percent, fill=Con)) +
    geom_bar(stat='identity') +
    labs(x='Paralog Status',y='Percentage',title = paste('DE ',SETTER,' Conservation Breakdown By Paralog Status', sep='')) +
    scale_fill_manual(values = Plate) +
    coord_flip() +
    scale_x_discrete(limits=c('NeoF','Other')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom") +
    geom_text(aes(label = sprintf("%1.2f%%", percent), y = position))
  print(p)
  
  p <- ggplot(G11, aes(x=joint, y=percent, fill=Con)) +
    geom_bar(stat='identity') +
    labs(x='Peak Classification',y='Percentage',title = paste('NeoF-DE ',SETTER,' Conservation Breakdown By Peak Type', sep='')) +
    scale_fill_manual(values = Plate) +
    coord_flip() +
    scale_x_discrete(limits=orders1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom") +
    geom_text(aes(label = sprintf("%1.2f%%", percent), y = position))
  print(p)
  
  dev.off()
}  
  