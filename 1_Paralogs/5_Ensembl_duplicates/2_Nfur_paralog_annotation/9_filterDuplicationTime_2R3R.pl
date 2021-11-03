# Filter for 3R and 2R nodes

use strict;
use warnings;

# Outfile for Ohnolog runs
open OUT3R, '>8_final_nfur_duplicates/nfurzeri_ReconsiledParalogs_Filtered_3R.txt' or die $!; # 3R
print OUT3R "Id1	Id2\n";
open OUT2R, '>8_final_nfur_duplicates/nfurzeri_ReconsiledParalogs_Filtered_2R.txt' or die $!; # 3R
print OUT2R "Id1	Id2\n";

# Open combined file 
open FH, '8_final_nfur_duplicates/Final_nfurzeri_paralogs.txt' or die $!;
my @paralogs = <FH>;

foreach (@paralogs){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;

	# If the reconciled node is vertebrate print it
	if ($line[2] eq 'Vertebrata' || $line[2] eq 'Euteleostomi' || $line[2] eq 'Chordata' || $line[2] eq 'Sarcopterygii' || $line[2] eq 'Neopterygii'){
		print OUT2R "$line[0]\t$line[1]\n";
	}
	if ($line[2] eq 'Teleosts' || $line[2] eq 'Clupeocephala' || $line[2] eq 'Acanthomorphata'){
		print OUT3R "$line[0]\t$line[1]\n";
	}
		
}

print `date`;


