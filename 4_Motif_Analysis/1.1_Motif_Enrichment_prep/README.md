README 
	Pipeline for 
	Last edited on 10/13/2021 -- GAREEVES

PIPELINE OVERVIEW
	This pipeline takes in various list of genes and peaks of interest and produces output genomic sequences/coordinates for those sites across 
	several fish species. These coordinate lists can then be handed to HOMER to evaluate the enrichment of TF binding sites within each set.

HARDWARE SPECIFICATIONS	
	The pipeline is designed to work on a MACOSX device (minimum version BigSir 11.4)	

TEST CASE
	No test data provided (Can be generated using other pipelines)
	

------------------------------------------------------------------------------------------	
DIRECTORY SETUP									|
	./HOME/										| -- The directory for the pipeline
		/Data/									| -- The directory for input and output data
			/Dev_<LIST_TYPE>_<OVERLAP>/			| -- The directories developmental peak output
			/Dia_Relaxed_<OVERLAP>/				| -- A directory for additional developmental peak output
			/LowCV_<DATA_TYPE>_<OVERLAP>/		| -- The directory for low coverage genes and associated peak output
			/LowCV_GENES_DE_<OVERLAP>/			| -- The directory for low coverage genes and associated differential peak output
			/NEO_1_and_2/						| -- The directory for peaks associated with neo-functionalized paralogs
			/Singletons_<GROUP>_<OVERLAP>/		| -- The directories for singleton genes and associated peak output
			/Singletons_<GROUP>_DE_<OVERLAP/	| -- The directories for singleton genes and associated differential peak output
												|
REQUIRED FILES									|
	INPUT DATA									|
		<SPECIES>_annotation.txt				| - Anottation file for peaks from a give species, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
		<SPECIES>_rpkm.txt						| - Files with RPKM-normalized read counts for peaks, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
		Peak_Master_<OVERLAP>.txt				| - Files showing conservation of peaks acrosses species at various overlap rates, generated in ~/HOME/ATACseq_Analysis/2_Conservation/
		Singleton_<GROUP>.txt					| - Singleton gene list associated with various RNAseq states, generated in ~/HOME/Paralog_Classification/
		lowCV_GENES.txt							| - Genes of low coverage in RNAseq data, generated in ~/HOME/RNAseq_Analysis/
		Dev_<TYPE>.bed							| - Peaks differentially accessible during development, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
		Dia_Relaxed.bed							| - Peaks differentially accessible during diapause, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
		fish4_WGA.maf							| - Multiple, whole genome alignment between fish species, generated in ~/HOME/Multi_Alignment_Construction
		lowCV_PEAKS.bed							| - Peaks associated with low coverage genes from RNAseq analysis, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
		Paralog_List.txt						| - List of Paralog pairs in the African turquoise killifish and their expression state as various timepoints, generated in ~/HOME/Paralog_Classification/
		Nfur_Up_Master.bed						| - List of all peaks with differential chromatin accessibility during diapause, generated in ~/HOME/ATACseq_Analysis/2_Peak_Calling/
												|
	SCRIPTS										|
		NeoF_DE.py								| -> placed in ~/HOME/
		Peak_Cons_Bed_Generator_0.py			| -> placed in ~/HOME/
		Peak_Cons_Bed_Generator_1.py			| -> placed in ~/HOME/
		Peak_Cons_Bed_Generator_2.py			| -> placed in ~/HOME/
												|
												|
	SCRIPT RUN ORDER							|
												|
		NeoF_DE.py								|
												|
		Peak_Cons_Bed_Generator_0.py			|
		|										|
		V										|
		Peak_Cons_Bed_Generator_1.py			|
		|										|
		V										|
		Peak_Cons_Bed_Generator_2.py			|  
												|	
PRE-RUNNING INSTRUCTIONS						|
	Script updates								|
		NeoF_DE.py								| -> Path update
												|	Update the path to home directory on line 11
												|
		Peak_Cons_Bed_Generator_0.py			| -> Path update
												|	Update the path to home directory on line 12
												| -> Data list update
												|	If including extra data files add them to the list on line 13
												| 
		Peak_Cons_Bed_Generator_1.py			| -> Path update
												|	Update the path to home directory on line 12
												| -> Data list update
												|	If including extra data files add them to the list on line 13
												| -> Overlap list update
												|	If including extra peak overlap files add them to the list on line 14
												| 
		Peak_Cons_Bed_Generator_2.py			| -> Path update
												|	Update the path to home directory on line 14
												| -> Grouping list update
												|	If including extra group comparisons add them to the list on line 15
												| -> Data list update
												|	If including extra data files add them to the list on line 16
												|										
												|
												|
------------------------------------------------------------------------------------------
		