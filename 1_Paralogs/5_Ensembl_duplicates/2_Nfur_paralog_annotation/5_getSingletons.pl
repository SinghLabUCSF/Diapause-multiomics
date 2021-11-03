# get the genes that are in single copy in other fish i.e. don't have any paralog
# AND have a BBH wih nfurzeri. ---> I WILL CLASSIFY THEM AS SINGLETONES IN NFUR
# In the BBH file - the first column is for nfurzeri
# Change $other_fish_gene_col to 1 for drerio. For other fish it must be 0
#
use strict;
use warnings;

# Some variables -------------------------------
my $fish_to_compare = 'drerio';
my $nfur_gene_col = 0; # the colum number for the gene id for nfur. Since in the All genes file there are gene symbols for Nfurzer, I will take symbols here
my $other_fish_gene_col = 1; # Column number for the fish being compared. 
# 1 for drerio, 0 for everyone else


# ****** For drerio it's 1 for other fish it would be 0 **********

# Open output file -----------------------------
open OUT, ">6_Sigletons\/nfurzeri_singletons_from_$fish_to_compare\.txt" or die $!;
print OUT "nfurzeri singleton	$fish_to_compare BBH\n";

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

# Get the duplication time from paralog file for other fish -----------------
# Read the BBH file for the other fish and get the corresponding BBH in nfur
# if both the paralogs ahve a BBH, get the duplication time and print in outfile

my %Duplicate; # Ensembl paralogs are duplicates printed in both the directions, so duplicates can be avoided using this hash
foreach (</Volumes/Mybook_2/Ohnologs/Synteny_All_2016_03_09/2_Paralogs/1_Paralogs_from_7_latest_versions/$fish_to_compare\_gene_ensembl_*.archive.ensembl.txt>){
	
	print "$_\n";
	
	open PARALOG, "$_" or die $!;
	# ****** Comment lines related to Duplicate hash if all paralogs should be printed just like for Ensembl

	foreach (<PARALOG>){ # foreach pair
	
		my @line = split "\t", $_;
		map {$_ =~s/^\s+|\s+$|^\t|\t$|\n//g} @line;
	
		if ($line[1] eq 'NA' || $line[1] eq ''){ # If there is no paralog
		
			if (exists $BBH{$line[0]}){ # and the gene has a BBH
			
				if (not exists $Duplicate{$line[0]}){ # and it has not been checked before				

					$Duplicate{$line[0]} = ''; 				# make a hash
					print OUT "$BBH{$line[0]}\t$line[0]\n";
				
					#print OUT "$BBH{$line[0]}\t$BBH{$line[1]}\t$line[2]\t$line[0]\t$line[1]\n";
				}
			}
		}
	}
}

print `date`;



