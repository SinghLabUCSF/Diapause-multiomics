## RNA-seq data processing and classification of paralogs for Nothobranchius furzeri

1. [1_getCountsAndTPM](1_getCountsAndTPM.R)  
   To generate count matrix for each library.
   
2. [2_getCombinedCountMatrix](2_getCombinedCountMatrix.pl)  
   Combine all count matrices into a single file.
   
3. [3_DESeq2_nfur_diapause_CK.R](3_DESeq2_nfur_diapause_CK.R)   
   Do pairwise differential expression analysis.

4. [4_join_DE-lists_DESeq2](4_join_DE-lists_DESeq2.R)
   Combine all pairwise lists from step above and generate two lists of genes up and down in diapause.
   
5. [7b_classifyAllDuplicatesFC](7b_classifyAllDuplicatesFC.pl)   
   Classify all into categories and identify neofunctionalized paralogs.
   
6. [7c_reordeClassifiedDuplicates_Ensembl](7c_reordeClassifiedDuplicates_Ensembl.pl)   
   [7c_reordeClassifiedDuplicates_Orthofinder](7c_reordeClassifiedDuplicates_Orthofinder.pl)   
   Reorder paralogs based on their upregulation in diapause. For neofunctionalized paralogs, the gene listed first in the paralog pair is always upregulated in diapause.

7. [7d_addAge_Orthofinder-71](7d_addAge_Orthofinder-71.pl)   
   [7d_addAge_Ensembl.](7d_addAge_Ensembl.pl)  
   [7d_addAge_Orthofinder-13](7d_addAge_Orthofinder-13.pl)  
   [7d_addAge_Orthofinder-31](7d_addAge_Orthofinder-31.pl)   
   Add the age category or duplication time of paralog to the file.
   
8. [8_Node_counts_71.pl](8_Node_counts_71.pl)
   [8_Node_counts_31.pl](8_Node_counts_31.pl)
   [8_Node_counts_13.pl](8_Node_counts_13.pl)   
   [8_Node_counts_Ensembl.pl](8_Node_counts_Ensembl.pl) 
   [8_Node_counts_ancientSSD.pl](8_Node_counts_ancientSSD.pl)  
   To counts the paralogs beonging to each duplication time. Very ancient (vertebrates), ancient (fish), recent/very recent (killifish) that are specialized for diapause.   
   I am making aseparate script for each because the nide names can be different for each category.
   
9. [8_Node_counts_bootstrap_Orthofinder-71.pl](8_Node_counts_bootstrap_Orthofinder-71.pl)   
   [8_Node_counts_bootstrap_Orthofinder-31.pl](8_Node_counts_bootstrap_Orthofinder-31.pl)   
   [8_Node_counts_bootstrap_Orthofinder-13.pl](8_Node_counts_bootstrap_Orthofinder-13.pl)   
   [8_Node_counts_bootstrap_Ensembl.pl](8_Node_counts_bootstrap_Ensembl.pl)   
   [8_Node_counts_bootstrap_Orthofinder-71_ancientSSD.pl](8_Node_counts_bootstrap_Orthofinder-71_ancientSSD.pl)   
   To do 10000 bootstrap of the datasets and count the ancient/recent/very recent paralogs that are specialized for diapause. I will plot them in prism for the paper.   
   I am making aseparate script for each because the nide names can be different for each category.
   
   > The results of these are used to plot the specialized age-distributions throughput the paper.
   
10. [global_plots.R](global_plots.R) Plot the boxplot for expression of paralogs.

11. [ChiseqTest.R](ChiseqTest.R)   
    Do Chi-square tests for all the counts that are in the paper for N. furzeri and A. limnaeus.
   