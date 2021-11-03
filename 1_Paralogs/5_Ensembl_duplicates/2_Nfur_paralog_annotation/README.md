### Get the BBH between fish and N. furzeri

*Get the BBH to get duplication time for N. furzeri genes.*

1. Run blast between nfurzeri and 5 other fish in both the directions.
   **`run_BLASTp_Seq-database.pl`**
   All these are with e-value 1-003
   
2. I modify the file names to sensible names. Originals are in [1_original_files](2_best_hit_files) folder.
   All Best hit output files are in [2_best_hit_files](2_best_hit_files)

3. Identify BBH between nfureri and all other fish. For all the nfur proteins I print it's
   BBH, otherwise the column will be black.
   **`1_Identify_BBH.pl`**
   [3a_bbh_outfiles](3a_bbh_outfiles)
   [3b_all_hit_outfiles](3b_all_hit_outfiles)
   
4. **`2_getDuplicationTime.pl`** gets the duplication time for BBH genes for 5 different fish for 6 different versions.
   In [4_get_duplication_time_ensembl](4_get_duplication_time_ensembl)   

5. **`3_getAllNodes.pl`**    
   Reconcile the duplication node from different ensembl versions.  
   Outfile [5_Combined_nodes](5_Combined_nodes)
   
6. **`4_getConsensusDuplicationNodes.pl`** To get all duplication nodes in N. furzeri based on all fish. And then consensus
   duplication nodes.
   Output: In [5_Combined_nodes](5_Combined_nodes)   
   Final files: _nfurzeri_all_dup_nodes.txt_ and _nfurzeri_consensus_dup_nodes.txt_

7. **`5_getSingletons.pl`**    
    To get the singletons from 5 fish and all 6 Ensembl versions.
    Outfile: [6_Sigletons](6_Sigletons)

8. **`6_getParalogsFromFamilies.pl`**      
   Outfile: [7_get_duplication_time_from_families](6_get_duplication_time_from_families)   
   This is to get the paralogs that are in the gene family clustering but not in Ensembl. This will get the paralogs duplicated
   in the fish that are not part of Ensembl pipeline. 
   Output _nfur_duplicates_from_families.txt_ See _Nodes.JPG_ for the nodes manually combined from Ensembl and NCBI Taxonomy.

9. **`7_combineParalogs.pl`** and  **`8_combineSingletons.pl`**
    Combine paralogs and singletons from both gene families and Ensembl 
    Outfile: [8_final_nfur_duplicates](8_final_nfur_duplicates) 

10. **`9_filterDuplicationTime_2R3R.pl`** Filter duplication time to run OHNOLOG 
   _nfurzeri_ReconsiledParalogs_Filtered_2R/3R.txt_ in [8_final_nfur_duplicates](8_final_nfur_duplicates) directory.   
   Copy them for ohnologs run to: /Volumes/Mybook_2/Ohnologs/Synteny_All_2016_03_09/2_Paralogs/4_filtered_paralogs
   
