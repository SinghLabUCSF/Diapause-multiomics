# To identify BBH between the fish genomes
# All the nfurzeri gene would not be there because all genes do not have a hit
# BBH will always have same number of gene in both the genomes

use strict;
use warnings;

# Run one by one for multiple genomes
my $from = 'nfurzeri';
my $to = 'olatipes';

open FORWARD, "2_best_hit_files\/BestHits_$from\-to-$to\.txt" or die $!;
open REVERSE, "2_best_hit_files\/BestHits_$to\-to-$from\.txt" or die $!;
open BBH, ">3a_bbh_outfiles\/BBH_$from\-$to.txt" or die $!;
open ALLHITS, ">3b_all_hit_outfiles\/All_Hits_$from\-$to.txt" or die $!;

# hash for key value pair
my %firstK_Vpair;

# construct the hash for one of the file - this is file for second column
foreach (<REVERSE>){
	
	my @line = split " ", $_;
	$line[0] =~s/^\s+|\s+$|^\t|\t$|\n//g;
	$line[1] =~s/^\s+|\s+$|^\t|\t$|\n//g;
	
	my $comp = $line[0].'__'.$line[1];		
	$firstK_Vpair{$comp} = \@line;
}

# For each pair in the second file - for the first column
foreach (<FORWARD>){
	
	my @line = split "\t", $_;
	$line[0] =~s/^\s+|\s+$|^\t|\t$|\n//g;
	$line[1] =~s/^\s+|\s+$|^\t|\t$|\n//g;
	
	my $comp = $line[1].'__'.$line[0];
	
	if (exists $firstK_Vpair{$comp}){		
		print BBH "$line[0]\t$line[1]\n";
		print ALLHITS "$line[0]\t$line[1]\tY\n";
	
	}
	else {
		print ALLHITS "$line[0]\t$line[1]\tN\n";
	}
}

print `date`;
