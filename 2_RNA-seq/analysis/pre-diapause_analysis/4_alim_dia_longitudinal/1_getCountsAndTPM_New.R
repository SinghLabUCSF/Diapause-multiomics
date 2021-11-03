# To do differential enrichment analysis by combining multiple samples
# Check what is the attribute type in your gtf file = "gene_name"
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

# Input variables ---------------------------
genename = "gene_name"
gtffile = "/Volumes/Mybook_3/Other_organism_genomes/alim/GFF/ref_Austrofundulus_limnaeus-1.0_top_level.gtf" # This needs to download from NCBI
outfilepath = "../../../feature_counts"
sampleTable <- read.csv("../../../fastq/alim_dia_longitudinal/experiment_design_alim.csv", row.names = 1, stringsAsFactors = F)
# ------------------------------------------
sampleTable$files

# Should be all TRUE
file.exists(sampleTable$files)
file.exists(gtffile)
file.exists(outfilepath)

library("Rsubread") # Use to get feature counts and gene lengths

for (i in 1:length(sampleTable$files)){
  
  outfilename <- file.path(outfilepath, paste0(sampleTable$lib[i], "_counts.csv"))
  print (sampleTable$files[i])
  outfilename
  cts = featureCounts(sampleTable$files[i], annot.ext = gtffile, isGTFAnnotationFile = T, GTF.attrType = genename)
  # Take data in a temporary variable, and change the column name to a meaningful one
  x = data.frame(cts$annotation[,c("GeneID","Length")],cts$counts,stringsAsFactors=FALSE)
  colnames(x)[3]
  colnames(x)[3] = "Counts"
  
  # get tpm, this is not normalized across samples
  x$TPM = tpm(x$Counts, x$Length)  
  write.table(x,file=outfilename,quote=FALSE,sep=",",row.names=FALSE)

  rm(outfilename, cts, x) # remove the temp variables for the next run

}

