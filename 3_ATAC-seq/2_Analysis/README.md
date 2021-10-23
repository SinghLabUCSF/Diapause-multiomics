## ATAC-seq data processing and analysis


1. [DiffBind_Deseq2.R](DiffBind_Deseq2.R)  
   [DiffBind_EdgeR.R](DiffBind_EdgeR.R)  
   These is the script to do differential abundance analysis, for all pairwise contrasts, for all conditions, for all species.  
   This was done using DiffBind version 2.16.2.
   It also generates a consensus peak file that I use for my cross species comparison. These results are saved in DeSeq2 and EdgeR folders for N. furzeri.
   
2. [nfur_DEseq2](getMasteFile_DEseq2.pl)    
   [nfur_EdgeR](getMasteFile_EDGER.pl)  
   To get a master peak file that includes all conditions and information.
   
3. Combine results from both DESeq2 and EDGER.   
   * [nfur_final_peaks/1_getCombinedMasterFile_All.pl](nfur_final_peaks/1_getCombinedMasterFile_All.pl)   
     Combine peaks from both EDGER and Deseq2.  
        
   * [4_getListsDiapause.pl](4_getListsDiapause.pl)  
     Script to select peaks that change between development and diapause but not between the two development conditions.   
     
   * [UpSetPlot.R](UpSetPlot.R)  Plot an Up Set plot for these peaks.
   
4. [annotate-peaks_ChIPseeker.R](annotate-peaks_ChIPseeker.R)   
   Annotate the peaks and plot annotaton as a bar plot.
   
5. Get correlation between RNA-seq and ATAC-seq
   * [correlate_with_RNA_overall_new.pl](correlate_with_RNA_overall_new.pl) Get median and maximum ATAC and RNA values per gene.   
   * [correlate_plot.R](correlate_plot.R) plot the correlation between medians.   
   