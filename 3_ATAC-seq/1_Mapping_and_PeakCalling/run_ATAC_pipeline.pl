use strict;
use warnings;

# USAGE NOTES: Chekc if the $fastqpath has just the names or full paths for input files, 
# and adjust the paths in command accordingly. E.g. for nfur, there are full paths and 
# for alim, only the file names.
# Also see what's the pattern for library name at line 33 onwards and change.
#

# -----------  Aaul development final run 1 ---------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul/files_aaul.txt';    # WITHOUT \n AT THE END. 
#my $in_dir = "fastq/nfur_aaul";                       # Input dir for fastq files to start with. Use full path
#my $out_dir = "aaul_run1";                            # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/aaul/aaul';            # BT2 index directory for this oragnism
#my $mitochr = 'REMOVED';                              # mito chromosome to be removed. Already removed for Aaul.       
#my $gene_coordinates = 'Aaul_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 838719370;                # for macs2: Nfur:856808511 Aaul:838719370
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Aaul and Ast development final run 2 ----------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_ast/aaul_ast.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_ast";              # Input dir for fastq files to start with. Use full path
#my $out_dir = "aaul_ast_run2";                        # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/aaul/aaul';            # BT2 index directory for this oragnism
#my $mitochr = 'REMOVED';                              # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Aaul_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 838719370;                # for macs2: Nfur:856808511 Aaul:838719370
#------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------   Aaul genome pilot 3 ----------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_pilot3/files_aaul.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_pilot3";            # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_aaul_alim_pilot3";                 # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/aaul/aaul';             # BT2 index directory for this oragnism
#my $mitochr = 'REMOVED';                               # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Aaul_gene_coordinates.bed';    # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 838719370;                 # for macs2: Nfur:856808511 Aaul:838719370 
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Alim development pilot 3 ----------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_pilot3/alim_dev.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_pilot3";            # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_aaul_alim_pilot3";                 # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/alim/alim';             # BT2 index directory for this oragnism
#my $mitochr = 'REMOVED';                               # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Alim_gene_coordinates.bed';    # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 483220538;                 # for macs2: Nfur:856808511 Aaul:838719370 
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Nfur genome pilot 3 ----------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_pilot3/nfur_1m_dia.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_pilot3";            # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_aaul_alim_pilot3";                 # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/nfur/nfur';             # BT2 index directory for this oragnism
#my $mitochr = 'NC_011814.1';                           # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Nfur_gene_coordinates.bed';    # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 856808511;                 # for macs2: Nfur:856808511 Aaul:838719370 
# ------------------------------------------------------------------------------------------------------------------------------------------------------


# -----------  Nfur final run 1 ----------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul/files_nfur.txt';    # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul";                       # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_run1";                            # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/nfur/nfur';            # BT2 index directory for this oragnism
#my $mitochr = 'NC_011814.1';                          # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Nfur_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 856808511;                # for macs2: Nfur:856808511 Aaul:838719370
#-----------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Nfur final run 2 ----------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_ast/nfur.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_ast";              # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_run2";                            # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/nfur/nfur';            # BT2 index directory for this oragnism
#my $mitochr = 'NC_011814.1';                          # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Nfur_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 856808511;                # for macs2: Nfur:856808511 Aaul:838719370
#-----------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Alim final run 2 ----------------------------------------------------------------------------------------------------------------------
my $fastq_list = 'fastq/nfur_aaul_alim_ast/alim.txt';    # WITHOUT \n AT THE END. 
my $in_dir = "fastq/nfur_aaul_alim_ast";                       # Input dir for fastq files to start with. Use full path
my $out_dir = "alim_run2";                            # outfile for fastqc, fastq, read alignment and tpm
my $genome_dir = 'genome_index/alim/alim';            # BT2 index directory for this oragnism
my $mitochr = 'UNKNOWN';                              # mito chromosome to be removed. Already removed for Aaul.       
my $gene_coordinates = 'Alim_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
my $effective_genome_size = 483220538;                # for macs2: Nfur:856808511 Aaul:838719370 Alim:483220538
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Drer from nature paper  ----------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/Drer/zebrafish_files.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/Drer";            # Input dir for fastq files to start with. Use full path
#my $out_dir = "drer_nature";                 # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/drer/drer';             # BT2 index directory for this oragnism
#my $mitochr = 'MT';                           # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Drer_gene_coordinates.bed';    # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 1368780147;                 # for macs2: Nfur:856808511 Aaul:838719370 
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------  Olati from nature paper  ----------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/Olat/medaka_files.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/Olat";            # Input dir for fastq files to start with. Use full path
#my $out_dir = "olat_nature";                 # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/olat/olat';             # BT2 index directory for this oragnism
#my $mitochr = 'MT';                           # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Olat_gene_coordinates.bed';    # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 733566086 ;                 # for macs2: Nfur:856808511 Aaul:838719370 
# -----------------------------------------------------------------------------------------------------------------------------------------------------

# ------------ All librries aligned to Nfur ------------------------------------------------------------------------------------------------------
#my $fastq_list = 'fastq/nfur_aaul_alim_ast/nfur.txt'; # WITHOUT \n AT THE END.
#my $in_dir = "fastq/nfur_aaul_alim_ast";              # Input dir for fastq files to start with. Use full path
#my $out_dir = "nfur_run2";                            # outfile for fastqc, fastq, read alignment and tpm
#my $genome_dir = 'genome_index/nfur/nfur';            # BT2 index directory for this oragnism
#my $mitochr = 'NC_011814.1';                          # mito chromosome to be removed. Already removed for Aaul.
#my $gene_coordinates = 'Nfur_gene_coordinates.bed';   # Gene coordinates for TSS heatmap: Nfur_gene_coordinates.bed, Aaul_gene_coordinates.bed
#my $effective_genome_size = 856808511;                # for macs2: Nfur:856808511 Aaul:838719370
#-----------------------------------------------------------------------------------------------------------------------------------------------------


#if (!-e "fastq/$out_dir"){print `mkdir -p fastq/$out_dir`}
if (!-e "fastqc/$out_dir"){print `mkdir -p fastqc/$out_dir`}
if (!-e "alignment_SAM/$out_dir"){print `mkdir -p alignment_SAM/$out_dir`}
if (!-e "alignment_BAM/$out_dir"){print `mkdir -p alignment_BAM/$out_dir`}
if (!-e "TSS_heatmaps/$out_dir"){print `mkdir -p TSS_heatmaps/$out_dir`}
if (!-e "peaks/$out_dir"){print `mkdir -p peaks/$out_dir`}
if (!-e "QC/$out_dir"){print `mkdir -p QC/$out_dir`}

open FH, "$fastq_list" or die $!; 

foreach (<FH>){
	
	my @line = split " ", $_;
	map {$_=~s/\n//g} @line;

	#$line[0] =~/(.+)_S\d.+\.fastq\.gz/g;           # for my nfur diapause data run1
	#$line[0] =~/(.+_S\d).+\.fastq\.gz/g;           # for my nfur diapause data run2
	#$line[0] =~/(.+_S\d+).+\.fastq\.gz/g;          # Aaul run 1
	$line[0] =~/(.+_S\d+)_.+\.fastq\.gz/g;          # Alim 
	#$line[0] =~/(.+)_\d+.fastq.gz/g;i              # zebrafish medaka
	
	my $libname = $1;
	print "Processing $libname ... \n";
	#                                 $1        $2       $3       $4          $5       $6       $7       $8                $9
 	#print `./ATAC_template.sh     $out_dir $line[0] $line[1] $genome_dir $libname  $in_dir  $mitochr $gene_coordinates $effective_genome_size`; # test without job submission
	#print "sbatch ATAC_template.sh $out_dir $line[0] $line[1] $genome_dir $libname $in_dir $mitochr $gene_coordinates $effective_genome_size\n";	
 	print `sbatch ATAC_template.sh $out_dir $line[0] $line[1] $genome_dir $libname  $in_dir $mitochr $gene_coordinates $effective_genome_size`;
 	#exit;
}

