README 
	Pipeline for comparing ATACseq and multi-species whole-genome alignment
	Last edited on 09/28/2021 -- GAREEVES

PIPELINE OVERVIEW
	This pipeline is used to generate several outputs and compile several metadata files comparing data between
	the whole-genome multi-alignment (~/Multi_Alignment_Construction/) and multiple species ATAC-seq data 
	(~/ATACseq_Analysis/1_Mapping/ & ~/ATACseq_Analysis/2_Peak_Calling/). Specifically the coordinates of 
	ATAC-seq peaks in various species are mapped to the blocks of the multi-alignment and then converted to the
	coordinates of the reference species (Nfur). This data is then used to make master lists showing corresponding
	peaks mapped to the same coordinates at various stringency (single-base overlap, 25% overlap, & 50% overlap) and
	do the same for each peak sets feature types (Intergenic, UTR, Promoter, Intron, Exon), Expression. These are
	then also used to generate metadata tables that contain peak info, feature info, paralog status 
	(Diapause-Specialized, Development-Specialized, Non-Specialized, Singletons), Paralog age (Very Ancient, Ancient,
	Recent), and closest genes. General stats are also calculated for the mutli-alignment including block-size breakdown
	and genome coverage in the alignment.
	
HARDWARE SPECIFICATIONS	
	The pipeline is designed to work on a MACOSX device (minimum version BigSir 11.4)

TEST CASE
	No test data provided (Can be generated using other pipelines)

------------------------------------------------------------------------------------------	
DIRECTORY SETUP						|
	./HOME/							| -- The directory for the pipeline
		/Data/						| -- The directory for input and output data
		/Packages/					| -- For required utilities (may not contain subdirectories, see below)
									|
									|
REQUIRED PACKAGES					|
	bedtools						| -- https://github.com/arq5x/bedtools2 (does not allow subdirectories)
		bedtools sort				|
		bedtools intersect			|
									|
REQUIRED FILES						|
	INPUT DATA						|
		<Species>_rpkm.txt			| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		<Species>_annotations.txt	| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		Nfur.<Species>.sing.maf		| -> Generated in ~/Multi_Alignment_Construction/, place in ~/HOME/Data/
		fish4_WGA.maf				| -> Generated in ~/Multi_Alignment_Construction/, place in ~/HOME/Data/
		Master_DE_Dia_Down_Peak.bed	| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		Master_DE_Dia_Down_Peak.bed	| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
									|
	SCRIPTS							|
		ATACseq_Conservation.py		| -> placed in ~/HOME/
									|	
PRE-RUNNING INSTRUCTIONS			|
	Path update						| if paths to all packages are in your .bashrc/.bash_profile,
									|	Search and remove all instances of '''PACK + ''' form ATACseq_Conservation.py
									| else:
									|	only follow the script updates below
	Script updates					|
		ATACseq_Conservation.py			| Update the current HOME directory path on line 22
							| if required packages are set to PATH do nothing, otherwise update path to packages on line 23
									| Adjust header names and data field count for 'RPKM_master' on lines 29-48
									| Adjust effective genomes sizes for on lines 54-58
									|
									|
------------------------------------------------------------------------------------------
		