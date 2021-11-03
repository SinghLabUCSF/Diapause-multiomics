# get the duplication time from other fish if both genes in paralog pairs have a BBH wih nfurzeri.
# In the BBH file - the first column is for nfurzeri
# Change $other_fish_gene_col to 1 for drerio. For other fish it must be 0
#
use strict;
use warnings;

# Some variables -------------------------------
my $fish_to_compare = 'drerio';
my $nfur_gene_col = 0; # the colum number for the gene id for nfur. Since in the All gnee file there are gene symbols for Nfurzer, I will take symbols here
my $other_fish_gene_col = 1; # Column number for the fish being compared. 
# ****** For drerio it's 1 for other fish it would be 0

# Read the BBH file ----------------------------
open FISH, "3a_bbh_outfiles\/BBH_nfurzeri-$fish_to_compare\.txt" or die $!;

my %BBH; # hash to store the BBH genes. Key => other fish; value => nfur BBH

foreach (<FISH>){
	
	my @line = split "\t", $_;
	map {$_ =~s/^\s+|\s+$|^\t|\t$|\n//g} @line;
	
	# Column 0 has nfur id, 1 has id for the other fish
	my @nfur = split '\|', $line[0]; 
	my @other_fish = split '\|', $line[1];
	
	my $nfurid = $nfur[$nfur_gene_col];
	my $other_fish_id = $other_fish[$other_fish_gene_col];
	
	#print "$nfurid\t$other_fish_id\n";
	$BBH{$other_fish_id} = $nfurid;
}
#print scalar keys %BBH;


foreach (</Volumes/Mybook_2/Ohnologs/Synteny_All_2016_03_09/2_Paralogs/1_Paralogs_from_7_latest_versions/$fish_to_compare\_gene_ensembl_*.archive.ensembl.txt>){
	
	print "$_\n";
	
	$_ =~/.+$fish_to_compare\_gene_ensembl_(.+)\.archive\.ensembl.txt/g;
	my $version = $1;
	print "$version\n";

	# Open output file -----------------------------
	open OUT, ">4_get_duplication_time_ensembl\/nfurzeri-paralogs_$fish_to_compare\_$version\.txt" or die $!;
	print OUT "nfurzeri paralog 1	nfurzeri paralog 2	duplication time	$fish_to_compare BBH for 1	$fish_to_compare BBH for 2\n";


	# Get the duplication time from paralog file for other fish -----------------
	# Read the BBH file for the other fish and get the corresponding BBH in nfur
	# if both the paralogs have a BBH, get the duplication time and print in outfile. Do this for all the fish species

	open PARALOG, "$_" or die $!;
	my %Duplicate; # Ensembl paralogs are duplicates printed in both the directions, so duplicates can be avoided using this hash
	# ****** Comment lines related to Duplicate hash if all paralogs should be printed just like for Ensembl

	foreach (<PARALOG>){ # foreach pair
	
		my @line = split "\t", $_;
		map {$_ =~s/^\s+|\s+$|^\t|\t$|\n//g} @line;
	
		if ($line[1] ne 'NA' || $line[1] ne '' || $line[2] ne 'NA' || $line[2] ne ''){ # If the pair has both the genes
		
			if ((exists $BBH{$line[0]}) && (exists $BBH{$line[1]})){ # and both are BBH
			
				if (not exists $Duplicate{$line[0]}{$line[1]}){ # and do not exists in Duplicate hash
				
					# make a hahs in both the directions
				
					$Duplicate{$line[0]}{$line[1]} = '';
					$Duplicate{$line[1]}{$line[0]} = '';
				
					print OUT "$BBH{$line[0]}\t$BBH{$line[1]}\t$line[2]\t$line[0]\t$line[1]\n";
				}
			}
		}
	}
	close (OUT);
	close (PARALOG);
}


print `date`;

