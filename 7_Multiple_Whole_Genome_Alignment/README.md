## Construction of multiple whole genome alignment

Pipeline for the construction of a multi-species whole-genome alignment

### PIPELINE OVERVIEW
This pipeline takes a set of genome fasta files and generates both pairwise alignment	and multi-alignment files. The genomes are split into individual chromosomes and then	pairwise alignments are generated in an all-by-all manner for the selected species.	These pairwise alignments are then used to build chromosome wide chains that are 	merged into a single, genome-wide chain and net. The pipeline then can use the final	pairwise alignments to generate various multiple alignments: (a) generate a reference	base, one-way multi-alignment, (b)generate a reference free, all-way multi-alignment,	or (c) generate a reference based, all-way multi-alignment. The pipeline can also 	clean up any intermediate files generated during the multi-alignment process.

### HARDWARE SPECIFICATIONS	
The pipeline is designed to work on a cluster using SLURM for job submission and a	86x64 ubuntu/linux architecture.

### TEST CASE
For the purpose of testing the pipeline a set of small test genomes (Exam, Test, and	Samp) are already included in ~/Data/1_genomes/Whole/ as well as an example formate	file at ~/Sub_Scripts/Species_List.txt. This example case should have and approximate	runtime of 30 minutes. 

------------------------------------------------------------------------------------------	
### DIRECTORY SETUP		 
	./HOME/					| -- The directory for the pipeline
		/Data/				| -- The directory for input and output data
			/1_genomes/		| -- For input genomes files
				/sizes/		| -- In which CHR sizes are calculated/stored
				/split/		| -- Into which fasta files are split
				/w2b/		| -- In which fasta files are stored as 2bit references
				/whole/		| -- In which the user places genome input fasta files
			/2_2bit/		| -- For 2bit chromosome files
			/3_lav/			| -- For pairwise alignment files
				/err/		| -- For alignment job stderr
				/out/		| -- For alignment job stdout
			/4_psl/			| -- For psl alignment files
			/5_chain/		| -- For chromosome alignment chains
				/err/		| -- For chaining job stderr
				/out/		| -- For chaining job stdout
			/6mchain/		| -- For genome alignment chains
			/7_pnet/		| -- For genome alignement prenets
			/8_net/			| -- For genome alignment nets
			/9_axt/			| -- For chained and netted genomes
			/10_maf/		| -- For final pairwise alignments
			/11_Final/		| -- For final multi-genome alignments
		/Packages/			| -- For required utilities (may not contain subdirectories, see below)
		/Sub_Scripts/		| -- For non-edited scripts that run pipeline
							
							
### REQUIRED PACKAGES			
	Lastz-distribution		| -- https://www.bx.psu.edu/miller_lab/
		lastz_D				|
	Multiz-TBA-distribution	| -- https://www.bx.psu.edu/miller_lab/
		multiz				|
		multic				|
		get_covered			|
		get_standard_headers|
		maf_checkThread		|
		maf_order			|
		maf_project			|
		maf_sort			|
		mafFind				|
		roast				|
		tba					|
	UCSC-Utils-distibution	| -- http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
		axtChain			|
		axtSort				|
		axtToMaf			|
		chainMergeSort		|
		chainNet			|
		chainPreNet			|
		faSize				|
		faSplit				|
		faToTwoBit			|
		lavToPsl			|
		netSyntenic			|
		netToAxt			|
							
### REQUIRED FILES		
	INPUT DATA				|
		genome1.fa			|
		...					| -> all genomes place in ~/HOME/Data/1_genomes/Whole/
		genomeN.fa			|
							|
	SCRIPTS					|
		Multi_Alignment.py	| -> placed in ~/HOME/
		Align_Submission.sh	| -> placed in ~/Home/Sub_Scripts/
		Chain_Submission.sh	| -> placed in ~/Home/Sub_Scripts/
		PSL_Submission.sh	| -> placed in ~/Home/Sub_Scripts/
							|
	SCORE MATRIX			|
		HoxD55.q			| -> placed in ~/HOME/Sub_Scripts/
							|
	Format File				|
		Species_List.txt	| -> placed in ~/HOME/
							| line1: comma seperated list of species (max 4 character str each)
							| line2: runtype followed by a '_' and species tree
							| Examples:
							|	GEN1,GEN2,...,GENn
							|	RUNTYPE_(GEN1 (GEN2 GENn))
							| RUNTYPE OPTIONS:
							| REFERENCE-BASED, ALL-WAY, BOTH 
							|	 (Ref-roast)    (tba) (Ref-tba)
							|[Currently only REFERENCE-BASED option functional]
							
							
### PRE-RUNNING INSTRUCTIONS
	Path update				| if paths to all packages are in your .bashrc/.bash_profile,
							|	Search and remove all instances of '''PACK + ''' form Multi_Alignment.py
							| else:
							|	follow script updates section below
	Script updates			|
		Multi_Alignment.py	| Update the current HOME directory path on line 16
							|
		Species_List.txt	| Add all needed species on line 1
							| Adjust runtype and tree on line 2
							|
		Alim.fa				|if Alim_old.fa exists, delete Alim.fa and rename Alim_old.fa to Alim.fa
							|
		Aaus.fa				|if Aaus_old.fa exists, delete Aaus.fa and rename Auas_old.fa to Aaus.fa
							
------------------------------------------------------------------------------------------
		