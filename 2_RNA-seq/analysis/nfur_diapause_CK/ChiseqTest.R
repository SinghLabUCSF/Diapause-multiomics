# This is to compute a table of Chi-square tests for all the conditions in the paper.
library(dplyr)

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Input file lists. For N fur and A.lim 
fileList = c( # Figure 1E, S3A
              "5_Node_counts/1_Orthofinder_71spp/All_Counts_NeoF_le20",
              "5_Node_counts/1_Orthofinder_71spp/All_Counts_NeoF_Only1DupEvent",
              "5_Node_counts/1_Orthofinder_71spp/All_Counts_NeoF_OnlyMostSimilar",
              "5_Node_counts/1_Orthofinder_71spp/All_Counts_NeoF_le20_AncientSSD",
              
              "5_Node_counts/2_Orthofinder_31spp/All_Counts_NeoF_le20",
              "5_Node_counts/3_Orthofinder_13spp/All_Counts_NeoF_le20",
              "5_Node_counts/4_Ensembl/All_Counts_NeoF_le20",
              
              "../alim_dia_nondia/NeoF_Counts_All_le20"
          )

# Check if all files are there, must be all TRUE
file.exists(paste0(fileList, ".txt"))

for (i in 1:length(fileList)){ # Foreach file
  
  values = read.table(file = paste0(fileList[i], ".txt"), sep = "\t", head = T) # Read the values
  
  # Do Chi-square test with reference to expected genome wide counts
  # An example of the values and command for Nfur vertebrate node is below.
  chiTest = values %>%
    rowwise() %>% 
    mutate(
      test_stat = chisq.test(c(NeoF_count,(NeoF_SUM - NeoF_count)),p=c(Ref_count/Ref_SUM,(Ref_SUM - Ref_count)/Ref_SUM))$statistic,
      p_val = chisq.test(c(NeoF_count,(NeoF_SUM - NeoF_count)),p=c(Ref_count/Ref_SUM,(Ref_SUM - Ref_count)/Ref_SUM))$p.value
    )
  chiTest = as.data.frame(chiTest) # Change data time from tibble to data frame
  row.names(chiTest) = row.names(values) # Reassign row names
  write.table(chiTest, file = paste0(fileList[i], "_withStats.txt"), quote = F, row.names = T, col.names = T, sep = "\t") # Write
}

# Here is an example test for 71 vertebrate very ancient category. P-values must match.
# chisq.test(c(5268,(6247 - 5268)),p=c(15737/20091,(20091 - 15737)/20091))$p.value

# This is for Unclassified because the header names are different. I should fix this to integrate better but don't have time.
unclassy = read.table(file = paste0("5_Node_counts/1_Orthofinder_71spp/All_Counts_Unclassified_le20", ".txt"), sep = "\t", head = T)

unclassychiTest = unclassy %>%
  rowwise() %>% 
  mutate(
    test_stat = chisq.test(c(Unclassified_count,(Unclassified_SUM - Unclassified_count)),p=c(Ref_count/Ref_SUM,(Ref_SUM - Ref_count)/Ref_SUM))$statistic,
    p_val = chisq.test(c(Unclassified_count,(Unclassified_SUM - Unclassified_count)),p=c(Ref_count/Ref_SUM,(Ref_SUM - Ref_count)/Ref_SUM))$p.value
  )
unclassychiTest = as.data.frame(unclassychiTest)
row.names(unclassychiTest) = row.names(unclassy)
write.table(unclassychiTest, file = paste0("5_Node_counts/1_Orthofinder_71spp/All_Counts_Unclassified_le20", "_withStats.txt"), quote = F, row.names = T, col.names = T, sep = "\t")


