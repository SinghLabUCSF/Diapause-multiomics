## This is to do GO enrichment analysis based on BBH list from zebrafish
#
# The starting object needs to be a data.frame with the 
# GO Id's in the 1st col, the evidence codes in the 2nd 
# column and the gene Id's in the 3rd
# 20151120: Changed to print another file with GO ids and genes
#

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("/Volumes/Mybook_2/Nfurzeri_PostGenome_Analysis/GO_KEGG_enrichment_NCBI")
#### IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####

library("GOstats")

# --------------------------------------------- INPUT FILES --------------------------------------------------------
frame = read.table(file ="1_GO_terms/GO_killifish-human-zebra_combined.txt", header = T, colClasses=c(rep("factor",3)))
# List of universe genes. Background.
#universe = read.table(file ="3_gene_sets/universe_all_genes.txt", header = T);
universe = read.table(file ="3_gene_sets/universe_ATACseq.txt", header = T);
# Genes to be tested
genes = read.table(file ="3_gene_sets/Relaxed_Any-Dia_peaks.txt", header = T);
# Output file
# ontology MF, BP, CC
ontolg = "BP"
outfilename = paste0("4_results/GO", ontolg, "_DiaUpAll_Peaks_NewBG.txt")
# Output file with genes in GO terms
gotermlist = paste0("4_results/Genes_GO", ontolg, "_DevelopmentUpAll_Peaks.txt") # IMPORTANT: Delete this file if it already exists
# Minimum number of genes for a term. Exclude anything less than these many genes.
mingenes = 5 # ********* THIS FILETR IS NOT BEING USED
# Maximum number of genes for a term to filter large general terms. Exclude anything more than this.
maxgenes = 1000 # ********* THIS FILETR IS NOT BEING USED
# Relative enrichment filter
relenrich = 0 # 0 will get everything
# ------------------------------------------------------------------------------------------------------------------

writeLines(paste0("Ontology: ", ontolg, "\n", "Enrichment outfile: ", outfilename, "\n", "Genes outfile: ", gotermlist, "\n", "Minimum genes: ", mingenes, "\n", "Maximum genes: ", maxgenes, "\n", "Relative enrichment: ", relenrich, "\n"))

#str(frame)

# This is just to get the 3 column. I already have these so I dont need it anyway!
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)

# put your data into a GOFrame object
goFrame=GOFrame(goframeData,organism="Human")
#head(goframeData)

# cast this object to a GOAllFrame object will tap into the GO.db package and populate this object with the implicated GO2All mappings for you
goAllFrame=GOAllFrame(goFrame)

# generate geneSetCollection objects
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# Process the universe list
universe = universe$id
universe = lapply(universe, as.character)
universe = unlist(universe)
head(universe)

# Process the gene list of interest
genes = genes$id
genes = lapply(genes, as.character)
genes = unlist(genes)
head(genes)

params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                             geneSetCollection=gsc, geneIds = genes, 
                             universeGeneIds = universe, 
                             ontology = ontolg,
                             pvalueCutoff = 1,
                             conditional = F, # To consider GO DAG structure or not 
                             testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over"

# call hyperGTest to do the test
Over <- hyperGTest(params)
head(summary(Over))

# calculate enrichment and add it to data frame.
# Relative enrichment factor (E-value) for a GO term = (count/size)/(size/universe)
enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))

# create a new frame
SummaryOver = data.frame(summary(Over), enrichment)
head(SummaryOver)

# I was trying a bunch of ways to reduce the number of GO terms. None of them work well so I am not doing them now
# Filter the Over variable on parameters other than P-value
# Filter the summary of OVER with size of the term, at least 2 genes for a go term
# Filter the Over variable on parameters other than P-value, plot and see a few of them
# plot(SummaryOver$OddsRatio)
# plot(SummaryOver$enrichment)
# mean(SummaryOver$OddsRatio[is.finite(SummaryOver$OddsRatio)])
# FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$Count <= maxgenes & SummaryOver$enrichment >= relenrich & SummaryOver$OddsRatio > mean(SummaryOver$OddsRatio[is.finite(SummaryOver$OddsRatio)])),]

# This is what a typical filter would look like
#FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]

# For now I am not filtering anything 
FilteredSummaryOver = SummaryOver
head(FilteredSummaryOver)

# adjust p value for multile correction
padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")

# Add padj to the data frame
FinalSummaryOver = data.frame(FilteredSummaryOver, padj)
head(FinalSummaryOver)

# write to a file
write.table(FinalSummaryOver, outfilename, quote = F, row.names = F, sep = "\t")

#rm(list = ls())

# --------------------- To get the genes for each Go terms ---------------------------------------------

# isolate indexes for the go terms in final results
#ind.GO <- is.element(names(Over@goDag@nodeData@data), FinalSummaryOver$GOBPID)
ind.GO <- is.element(names(Over@goDag@nodeData@data), eval(parse(text=paste("FinalSummaryOver$", "GO",ontolg,"ID", sep=''))))
selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]

# get a go terms and genes in a new vaiable for all the terms in the results of enrichment
goTerms <- lapply(selected.GO, 
                  function(x) x$geneIds)
names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]

# This will carete a new file "genesForGOTerms.txt" that will have GO terms and genes in each gO terms
# Genes can be duplicate, GO terms should not be
# Number of Go terms or lines should be equal to the enriched go terms as in the other file generated by this script
# This needs to be processed to generate the desired files

for (i in 1:length(goTerms)){
  
  test = as.data.frame(do.call(rbind, goTerms[i]))
  write.table(test, file = gotermlist, quote = F, col.names = F, append = T) # append each line - so make sure the file is empty for each run, or renamed after each run
  rm(test)
}

rm(list = (setdiff(ls(), c("goframeData", "goFrame", "goAllFrame", "gsc", "frame", "universe", "genes"))))
#rm(list = ls())

