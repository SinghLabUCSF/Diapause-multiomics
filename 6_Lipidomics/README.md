## Process Lipidomics data

1. [1_get_dominant_ion.pl](1_get_dominant_ion.pl)  
   To get ion counts for each lipid class. This will be used to select dominant ions for each class to filter later.   
   I filter this manually and geenrate a major ion file [embryo_major_ions.txt](embryo_major_ions.txt)
   
2. [2_filter.pl](2_filter.pl)  
   To filtere lipids based on major ions, internal standards or not annotated properly etc.  
   
3. [3_filter_duplicates_rows.pl](3_filter_duplicates_rows.pl)   
   Filter multiple lipids in different rows by taking the one with highest intensity.
   
4. [4_filter_duplicates_FA.pl](4_filter_duplicates_FA.pl)   
   For lipids with multiple fatty acids, filter them based on m-score and t-score.
   
5. [5_filter_duplicates_rows.pl](5_filter_duplicates_rows.pl)  
   Filtering multiple lipids can generate duplicate rows. Remove them using the same code as step above.
   
6. [6_normalize_wrt_standards.pl](6_normalize_wrt_standards.pl)  
   Normalize lipids with respect to internal standards.   
   
7. [7b_median_normalization_Phospholipids.R](7b_median_normalization_Phospholipids.R)  
   Normalize lipids across samples based on phospholipid concentrations.
   
8. [8_sum_class.pl](8_sum_class.pl)   
   Sum the class concentration to do class specific analysis, and also look at MUFA, PUFA and short/mediam/long/very-long-chain FA etc.

9. [8b_plot_correlation.R](8b_plot_correlation.R)   
   To plot correlation between phospholipids and proteins.
   
10. [9_plotPCA_NfurSamples.R](9_plotPCA_NfurSamples.R)   
    [9_plotPCA_Nfur_Ast.R](9_plotPCA_Nfur_Ast.R)   
    Plot PCAs.
   
11. [10_get_DE_lipids.R](10_get_DE_lipids.R)   
    [10b_get_DE_lipid_classes.R](10b_get_DE_lipid_classes.R)   
    To get differential lipids between conditions, for individual or class specific analysis.
    