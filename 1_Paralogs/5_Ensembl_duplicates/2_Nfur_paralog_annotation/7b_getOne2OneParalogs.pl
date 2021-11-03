# Check if there are duplicate pairs in reverse direction
# And get the number of duplicates for each gene

use strict;
use warnings;

open OUT, '>8_final_nfur_duplicates/Final_nfurzeri_paralogs_new.txt' or die $!;

open FH, '/Volumes/Mybook_3/Ohnologs/Synteny_All_2016_03_09/OhnologsForNfur/2_Paralogs/2_Nfur_paralog_annotation/8_final_nfur_duplicates/Final_nfurzeri_paralogs.txt' or die $!;
my @file = <FH>;
shift @file;

my %Duplicates;
my %Dupcounts; my %Dupnodes;

foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ((exists $Duplicates{$line[0]}{$line[1]}) || (exists $Duplicates{$line[1]}{$line[0]})){
		
		die "pair $line[0] - $line[1] exists already\n";
	}
	else {
		$Dupcounts{$line[0]} += 1;
		$Dupcounts{$line[1]} += 1;

		push @{$Dupnodes{$line[0]}}, $line[2];
		push @{$Dupnodes{$line[1]}}, $line[2];
	}
	
}

#foreach (keys %Dupcounts){
	
#	print "$_\t$Dupcounts{$_}\t";
#	print join ',', @{$Dupnodes{$_}};
#	print "\n";
#}

print OUT "Id1	Id2	Node	Source	Fish support (Ensembl node)	Total paralogs for Id1	Total paralogs for Id2	Sum total paralog pairs\n";
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	print OUT join ("\t", @line), "\t";
	print OUT "$Dupcounts{$line[0]}\t$Dupcounts{$line[1]}\t";
	print OUT eval ($Dupcounts{$line[0]} + $Dupcounts{$line[1]}), "\n";
	
	#if ($Dupcounts{$line[0]} == 1 && $Dupcounts{$line[1]} == 1){
	#	
	#	print OUT "$_";
	#}
	
}

print `date`;

