# Get the list of peaks to be used for analysis
use strict;
use warnings;

# Open master file and get the desired conditions EDGE-R
open MF, "MasterPeakFile_DEseq2-EdgeR_Combined_All.txt" or die $!;
my @masterfile = <MF>;
my $head = shift @masterfile;

open OUT1, ">Strict_All-Dia.bed" or die $!;
open OUT2, ">Intermediate1_6d-Dia.bed" or die $!;
open OUT3, ">Intermediate_1m-6d-Dia.bed" or die $!;
open OUT4, ">Relaxed_Any-Dia.bed" or die $!; # This is what I will use for our analysis. Change between any diapause and does not change between the two development conditions.

foreach (@masterfile){
	
	my @line = split "\t", $_;
	map {$_=~s/\t//g} @line;
	
	#print "$line[21]\t$line[24]\t$line[33]\t$line[36]\t$line[42]\t${$DEseq2{$line[3]}}[42]\n";
	
	# Strict condition
	if ($line[21] eq 'Dia1m.up' && $line[24] eq 'Dia1m.up' && $line[33] eq 'Dia6D.up' && $line[36] eq 'Dia6D.up' && $line[42] eq ''){
		print OUT1 "$_"
	}
	
	# Intermediate	
	if (($line[33] eq 'Dia6D.up' && $line[36] eq 'Dia6D.up') && $line[42] eq ''){
		print OUT2 "$_"
	}

	# Relaxed 1
	if (($line[21] eq 'Dia1m.up' || $line[24] eq 'Dia1m.up') && ($line[33] eq 'Dia6D.up' || $line[36] eq 'Dia6D.up') && ($line[42] eq '')){
		print OUT3 "$_"
	}

	# Relaxed 2
	if (($line[21] eq 'Dia1m.up' || $line[24] eq 'Dia1m.up' || $line[33] eq 'Dia6D.up' || $line[36] eq 'Dia6D.up') && ($line[42] eq '')){
		print OUT4 "$_"
	}
	
}
print `date`;
