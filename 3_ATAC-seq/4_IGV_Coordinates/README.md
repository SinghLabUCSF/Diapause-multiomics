README 
	Pipeline for obtaining IGV base pair anchors and window coordinates across species
	Last edited on 09/28/2021 -- GAREEVES

PIPELINE OVERVIEW
	This pipeline is used to generate several outputs and compile several metadata files comparing data between
	the whole-genome multi-alignment (~/Multi_Alignment_Construction/) and multiple species ATAC-seq data 
	(~/ATACseq_Analysis/1_Mapping/ & ~/ATACseq_Analysis/2_Peak_Calling/). Specifically the coordinates for 1) the 
	center of each ATACseq peak, 2)the start of each gene promoter, and 3) the relationship of identified paralog
	pairs in the reference species, N.fur. These coordinates are then translated to other fish species using the 
	whole-genome multi-alignment blocks and each species coordinate set is expanded to encompass a 4KB window 
	centered on the above obtained coordinate anchors. These windows are then outputted in pairs generated between
	either comparable peaks between paralogs, a paralog peak and its partner's promoter, or the promotors of two
	paralogs.
	
HARDWARE SPECIFICATIONS	
	The pipeline is designed to work on a MACOSX device (minimum version BigSir 11.4)

TEST CASE
	No test data provided (Can be generated using other pipelines)

------------------------------------------------------------------------------------------	
DIRECTORY SETUP						|
	./HOME/							| -- The directory for the pipeline
		/Data/						| -- The directory for input and output data
									|
									|
REQUIRED PACKAGES					|
	none							|
									|
REQUIRED FILES						|
	INPUT DATA						|
		Nfur_RPKM.txt				| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		Nfur_Peak_All.bed			| -> Generated in ~/ATACseq_Analysis/3_Conservation/, place in ~/HOME/Data/
		Nfur.gtf					| -> Placed in ~/HOME/Data/, acquired from --- NCBI
		Nfur_gene.bed				| -> Generated using NFUR.gtf, placed in ~/HOME/Data/
		Promoter_List.txt			| -> Manually curated list, placed in~/HOME/Data/
		Partner_List.txt			| -> Manually curated list, placed in ~/HOME/Data/
		fish4_WGA.maf				| -> Generated in ~/Multi_Alignment_Construction/, place in ~/HOME/Data/
		Peak_Master_0.0.txt			| -> Generated in ~/ATACseq_Analysis/3_Conservation/, place in ~/HOME/Data/
		Manual_Cans	.txt			| -> Manually curated list, place in ~/HOME/Data/
									|
	SCRIPTS							|
		ATACseq_IGV_Coordinates.py	| -> placed in ~/HOME/
									|	
PRE-RUNNING INSTRUCTIONS			|
	Script updates					|
		ATACseq_IGV_Coordinates.py	| Update the current HOME directory path on line 22
									| Adjust maximum maf block number on line 21
									|
									|
------------------------------------------------------------------------------------------
		