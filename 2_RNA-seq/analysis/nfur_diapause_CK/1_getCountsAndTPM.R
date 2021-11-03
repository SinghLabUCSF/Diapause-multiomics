#  _ __  __ _ _ _ __ _ _ __  
# | '_ \/ _` | '_/ _` | '  \ 
# | .__/\__,_|_| \__,_|_|_|_|
# |_| 
# 
# Get feature counts using Rsubread package in the feature_count directory. This also generates a column for feature length that I may use later for cross species comparison.
# Also calculate TPM is you need.
#

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
gtffile = "../../Datasets/Nfur_genome/ref_Nfu_20140520_top_level_filtered.gtf"
sampleTable <- read.csv("../../fastq/nfur_diapause_CK/experiment_design.csv", row.names = 1, stringsAsFactors = F) # Sample BAM files
# ------------------------------------------
sampleTable$files

#filenames <- file.path("/Volumes/Mybook_3/RNA-seq/read_alignment_star/nfur/nfur_diapause_CK/", paste0(sampleTable$lib, "_UniquelyMapped.bam"))
#file.exists(filenames) # Should be all TRUE
#filenames

file.exists(sampleTable$files)

library("Rsubread") # Use to get feature counts and gene lengths

#for (i in 1:length(sampleTable$files)){
for (i in 1:1){ # Test for 1 library. Uncomment the above sentence if you want to run for all.
  
  outfilename <- file.path("../../feature_counts", paste0(sampleTable$lib[i], "_counts.csv"))
  print (sampleTable$files[i])
  outfilename
  cts = featureCounts(sampleTable$files[i], annot.ext = gtffile, isGTFAnnotationFile = T, GTF.attrType = genename, isPairedEnd = T) # In earlier version it used to detect PE automatically.
  # Take data in a temporary variable, and change the column name to a meaningful one
  x = data.frame(cts$annotation[,c("GeneID","Length")],cts$counts,stringsAsFactors=FALSE)
  colnames(x)[3]
  colnames(x)[3] = "Counts"
  
  # get tpm, this is not normalized across samples
  x$TPM = tpm(x$Counts, x$Length)  
  write.table(x,file=outfilename,quote=FALSE,sep=",",row.names=FALSE)

  rm(outfilename, cts, x) # remove the temp variables for the next run

}
