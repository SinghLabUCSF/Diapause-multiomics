library('GenomicFeatures')
library(GenomicAlignments)
library(ChIPseeker)

# Setwd to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############# INPUT FILES #################
# Nfur
sample_sheet = 'nfur_EdgeR/sample_sheet_3.csv'
output_dir = 'nfur_EdgeR/'
gff_file = './ref_Nfu_20140520_top_level.gff3'

# Drer
#sample_sheet = 'drer/sample_sheet_drer.csv'
#output_dir = 'drer/'
#gff_file = '/Volumes/Mybook_3/Other_organism_genomes/drerio/Danio_rerio.GRCz11.100.gff3'

# Drer selected stages
#sample_sheet = 'drer_selected-stages/sample_sheet_drer.csv'
#output_dir = 'drer_selected-stages/'
#gff_file = '/Volumes/Mybook_3/Other_organism_genomes/drerio/Danio_rerio.GRCz11.100.gff3'

# Olati
#sample_sheet = 'olat/sample_sheet_medaka.csv'
#output_dir = 'olat/'
#gff_file = '/Volumes/Mybook_3/Other_organism_genomes/olati/Oryzias_latipes.ASM223467v1.100.gff3'

# Olati selected stages
#sample_sheet = 'olat_selected-stages/sample_sheet_medaka.csv'
#output_dir = 'olat_selected-stages/'
#gff_file = '/Volumes/Mybook_3/Other_organism_genomes/olati/Oryzias_latipes.ASM223467v1.100.gff3'

# Alim
#sample_sheet = 'alim/sample_sheet_alim.csv'
#output_dir = 'alim/'
#gff_file = '/Volumes/Mybook_3/Other_organism_genomes/alim/GFF/ref_Austrofundulus_limnaeus-1.0_top_level.gff3'

# Aaul-Ast combined
#sample_sheet = 'aaul/sample_sheet_aaul-ast.csv'
#output_dir = 'aaul-ast/'
#gff_file = '//Volumes/Mybook_3/Nfur_genomes_DRV/aaustrale/AAU.v1.2.gff'

# Aaul ONLY
#sample_sheet = 'aaul/sample_sheet_aaul-ast.csv'
#output_dir = 'aaul/'
#gff_file = '//Volumes/Mybook_3/Nfur_genomes_DRV/aaustrale/AAU.v1.2.gff'

# Ast ONLY
#sample_sheet = 'ast/sample_sheet_aaul-ast.csv'
#output_dir = 'ast/'
#gff_file = '/Volumes/Mybook_3/Nfur_genomes_DRV/aaustrale/AAU.v1.2.gff'

###########################################

# Making a TxDb object: https://seandavi.github.io/ITR/transcriptdb.html
# GTF file is giving me an error so I am using original ncbi GFF3
txdb = makeTxDbFromGFF(gff_file)

# Save txdb so I don't have to make it every time
#library(AnnotationDbi)
#saveDb(txdb, 'txdb.nfurncbi')

# Load txdb from file
# library(AnnotationDbi)
# txdb = loadDb(file = 'txdb.nfurncbi')
# genes(txdb)

# Annotate all peaks and write their annotation in a new file
consensus = readPeakFile(paste0(output_dir, "consensusPeaks.bed"))
annoList = annotatePeak(consensus, TxDb=txdb)
peakAnnotations = annoList@anno
write.table(mcols(peakAnnotations), file = paste0(output_dir, "consensusPeaks_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

plotAnnoBar(annoList)

############################################ Only for Nfur ###########################################################################
# # Test plotting coverage
# peakDia <- readPeakFile("/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/nfur/nfur_bt_very-sensitive_unique/NfurM1_peaks.narrowPeak")
# peakEsc <- readPeakFile("/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/nfur/nfur_bt_very-sensitive_unique/NfurE1_peaks.narrowPeak")
# peakDia
# 
# # weightCol is the heading of the column that you want to plot
# covplot(peak, weightCol="V5")
# covplot(peakDia, weightCol="V10", chrs=c("NC_029649.1", "NC_029650.1"))
# covplot(peakEsc, weightCol="V10", chrs=c("NC_029649.1", "NC_029650.1"))
# 
# covplot(peakDia, chrs=c("NC_029649.1", "NC_029650.1", "NC_029651.1","NC_029652.1","NC_029653.1","NC_029654.1","NC_029655.1","NC_029656.1","NC_029657.1","NC_029658.1"))
# covplot(peakEsc, chrs=c("NC_029649.1", "NC_029650.1", "NC_029651.1","NC_029652.1","NC_029653.1","NC_029654.1","NC_029655.1","NC_029656.1","NC_029657.1","NC_029658.1"))


# This section is to annotate all peaks for all the files ------------------
# Read sample description sheet. See help for the dba function for column specification
# Column names are fixed so it automatically reads the peaks etc.
# samples <- read.csv(paste0(output_dir,sample_sheet))
# samples
# 
# # Plot annotation for all the libraries
# files = as.list(samples$Peaks)
# names(files) = row.names(samples)

# peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
#                        tssRegion=c(-2000, 2000), verbose=FALSE)
# pdf(paste0(output_dir, "FeatureAnnotations_AllPeaks.pdf"))
# plotAnnoBar(peakAnnoList)
# dev.off()

# Annotate all the peaks for DE peak files  ----------------------------
#### THESE ARE UP BASED ON THE FIRST CONDITION
# If there is an erro here it's most likely due to 1st row in bed file bein a header row. Just delete that and it should run
DE_peaks <- list(D6D_vs_PreD = "nfur_EdgeR/Dia6D-Kvyoung_DE-up_EdgeR.bed",
                 D1m_vs_PreD = "nfur_EdgeR/Kvyoung-Dia1m_DE-down_EdgeR.bed",
                 D6D_vs_ND = "nfur_EdgeR/Dia6D-Escape_DE-up_EdgeR.bed",
                 D1m_vs_ND = "nfur_EdgeR/Escape-Dia1m_DE-down_EdgeR.bed")


DEpeaksAnno <- lapply(DE_peaks, annotatePeak, TxDb=txdb,
                      tssRegion=c(-2000, 2000), verbose=FALSE)

pdf(paste0(output_dir, "FeatureAnnotations_DEPeaks_selected_conditions.pdf"))
plotAnnoBar(DEpeaksAnno)
dev.off()

# # ANNOTATION FOR DIAPAUSE UP DOWN FILES
# # escape down is diapause up
# dpup = paste0(output_dir, "Escape-Dia1m_DE-down_EdgeR.bed")
# anno_dpup = annotatePeak(dpup, TxDb=txdb)
# peakAnnotations_dpup = anno_dpup@anno
# write.table(mcols(peakAnnotations_dpup), file = paste0(output_dir, "Escape-Dia1m_DE-down_EdgeR_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# # escape up is diapause down
# dpdown = paste0(output_dir, "Escape-Dia1m_DE-up_EdgeR.bed")
# anno_dpdown = annotatePeak(dpdown, TxDb=txdb)
# peakAnnotations_dpdown = anno_dpdown@anno
# write.table(mcols(peakAnnotations_dpdown), file = paste0(output_dir, "Escape-Dia1m_DE-up_EdgeR_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# 
# # ANNOTATION FOR AAUL HB condition
# aahb = "/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/aaul_ast_run2/Aaus_HBY_1_S4_peaks.narrowPeak"
# anno_aahb = annotatePeak(aahb, TxDb=txdb)
# peakAnnotations_aahb = anno_aahb@anno
# write.table(peakAnnotations_aahb, file = paste0(output_dir, "Aaus_HBY_1_S4_peaks_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# 
# # ANNOTATION FOR AST HB
# asthb = "/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/aaul_ast_run2/Astr_HBY_1_S9_peaks.narrowPeak"
# anno_asthb = annotatePeak(asthb, TxDb=txdb)
# peakAnnotations_asthb = anno_asthb@anno
# write.table(peakAnnotations_asthb, file = paste0(output_dir, "Astr_HBY_1_S9_peaks_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# 
# # ANNOTATION FOR MEDAKA STAGE 25
# medaka25 = "/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/olat_nature/Olat_st251_peaks.narrowPeak"
# anno_medaka25 = annotatePeak(medaka25, TxDb=txdb)
# peakAnnotations_medaka25 = anno_medaka25@anno
# head(peakAnnotations_medaka25)
# write.table(peakAnnotations_medaka25, file = paste0(output_dir, "Olat_st251_peaks_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# 
# # Zebrafish 8 somites
# zebra8 = "/Volumes/Mybook_3/ATAC_Seq_Killifish/peaks/drer_nature/Drer_8som1_peaks.narrowPeak"
# anno_zebra8 = annotatePeak(zebra8, TxDb=txdb)
# peakAnnotations_zebra8 = anno_zebra8@anno
# head(peakAnnotations_zebra8)
# write.table(peakAnnotations_zebra8, file = paste0(output_dir, "Drer_8som1_peaks_annotations.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# 
