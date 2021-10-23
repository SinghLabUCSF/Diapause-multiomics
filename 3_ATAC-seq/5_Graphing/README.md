README 
	Pipeline for
	Last edited on 09/29/2021 -- GAREEVES

PIPELINE OVERVIEW
	Figure_3B.R
	This pipeline is used to apply VST normalization to integer-converted RPKM values across species ATAC-seq
	data. After this conversion the script generates a variety of PCA graphs based on various comparisons sets
	between the provided species. It also after performing integer-conversion outputs a new RPKM file with the
	updated integer values. The script also specifically prints the PC2 loading for the comparison between the
	African turquoise killifish and the South American Killifish (Diapause/Development PC) to be used later
	by the script PC_Loading_Bed.py.
	
	Figure_3C_S10.R
	This pipeline is used to intake peak conservation data across multiple species and associated metadata and
	subdivide the data meaningfully by various stats (Diapause/Development, Species, Paralog Status, etc.). 
	After subdivision the script then converts discrete counts of each division of the data to percentages 
	pertaining to comparisons between various subdivisions. These percentages are then output as a data table
	and graphed as various bar-charts of percentage breakdown. This function is performed for both conservation
	of ATAC-seq peaks and conservation ('alignability') of their underlying regions. 
	
	PC_Loading_Bed.py
	This pipeline is used to generate bed files for peaks identified as contributing to PC2 
	(Diapause/Development PC) in the comparison of African turquoise killifish and South American killifish.	
	
HARDWARE SPECIFICATIONS	
	The pipeline is designed to work on a MACOSX device (minimum version BigSir 11.4)

TEST CASE
	No test data provided (Can be generated using other pipelines)

------------------------------------------------------------------------------------------	
DIRECTORY SETUP								|
	./HOME/									| -- The directory for the pipeline
		/Data/								| -- The directory for input and output data
											|											|
REQUIRED PACKAGES							|
	ggplot2									| -> https://ggplot2.tidyverse.org/
	dplyr									| -> https://dplyr.tidyverse.org/
	tidyverse								| -> https://www.tidyverse.org/
	DEseq2									| -> https://bioconductor.org/packages/release/bioc/html/DESeq2.html
											|
REQUIRED FILES								|
	INPUT DATA								|
		Nfur_RPKM.txt						| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		<species>_annotations.txt			| -> Generated in ~/ATACseq_Analysis/2_Peak_Calling/, place in ~/HOME/Data/
		Final_Con_List_<type>_Edit_R.txt	| -> Generated in ~/ATACseq_Analysis/3_Conservation/, place in ~/HOME/Data/
		<type>_Master_0.0.txt				| -> Generated in ~/ATACseq_Analysis/3_Conservation/, place in ~/HOME/Data/
		Nfur_Alim_PC2.txt					| -> Generated in ~/ATACseq_Analysis/5_Graphing/, place in ~/HOME/Data/
											|
	SCRIPTS									|
		Figure_3B.R							| -> placed in ~/HOME/
		Figure_3C_S10.R						| -> placed in ~/HOME/
		PC_Loading_Bed.py					| -> placed in ~/HOME/
											|	
PRE-RUNNING INSTRUCTIONS					|
	Script updates							|
		Figure_3B.R							| Update the current HOME directory path on line 5
											| 
		Figure_3C_S10.R						| Update the current HOME directory path on line 5
											| 
		PC_Loading_Bed.py					| Update the current HOME directory path on line 7
											| 
											|
											|
------------------------------------------------------------------------------------------
		