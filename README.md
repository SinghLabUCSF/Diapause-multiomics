## Diapause multi-omics

**This is the repository having code and data analysis pipeline for multi-omics study on diapause evolution.**

> #### Reference   
  Evolution of diapause in the African turquoise killifish by remodeling ancient gene regulatory landscape.  
  Param Priya Singh*, G. Adam Reeves*, Kévin Contrepois, Mathew Ellenberger, Chi-Kuo Hu, Michael P. Snyder, and Anne Brunet.   
  Preprint: https://www.biorxiv.org/content/10.1101/2021.10.25.465616v1   
  
  More information can be found in the directories below:

### Paralog identification and dating
-----
* [1_Paralogs](1_Paralogs):  Pipeline to identify all the paralogs and their duplication time relative to other animals. 
-----
### RNA-seq data analysis
* [2_RNA-seq](2_RNA-seq): Gene expression analyses in development and diapause, and their integration with paralogs to identify the paralog pairs that are specialized for diapause or development.
-----
### ATAC-seq data analysis
* [3_ATAC-seq](3_ATAC-seq): Pipeline to process chromatin accessibility data and integrate it with RNA-seq to identify the accessible regions that are specialized for diapause.
-----

### Motif enrichment analysis
* [4_Motif_Analysis](4_Motif_Analysis): Enrichment, integration and analysis of transcription factor binding motifs at paralogs across species.
-----
### Functional enrichment analysis
* [5_Functional_Enrichment](5_Functional_Enrichment): Functional enrichment analysis using Gene Ontology for ATAC-seq and RNA-seq datasets.
-----
### Lipidomics data processing and analysis
* [5_Lipidomics](5_Lipidomics): Pipeline to process and analyze the lipidomics data in development and diapause and across species.
-----
### Multiple Whole Genome Alignment
* [7_Multiple_Whole_Genome_Alignment](7_Multiple_Whole_Genome_Alignment): Pipeline to process and analyze the lipidomics data in development and diapause and across species.

