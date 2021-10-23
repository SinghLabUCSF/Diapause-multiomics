## Motif enrichment analysis and plotting

1. [1_runMotifs.pl](1_runMotifs.pl)  
   To run all the HOMER motif enrichment analysis. All the genomes were locally installed using HOMER install genome utility.   
   
2. [2_processMotifs_all.pl](2_processMotifs_all.pl)   
   Process motifs and remove redundany using TomTom motif clustering and manual curation.   
   
3. [3.1_combine_motifs_allPara.R](3.1_combine_motifs_allPara.R)   
   Combine and plot motifs enriched in diapause specific peak at paralogs.   
      
4. [3.5_combine_motifs_Singleton_Paralogs.R](3.5_combine_motifs_Singleton_Paralogs.R)   
   Combine and contrast motifs enriched in diapause specific peak at paralogs vs singleton genes.   

5. [3.7_combine_motifs_PositiveSelection.R](3.7_combine_motifs_PositiveSelection.R)   
   Combine and contrast motifs enriched in all diapause specific peak vs the ones in peaks under positive selection.      
   
