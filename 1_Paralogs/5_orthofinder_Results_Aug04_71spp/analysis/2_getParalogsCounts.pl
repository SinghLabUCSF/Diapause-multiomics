# Check if there are duplicate pairs in reverse direction e.g. A-B or B-A, and make a unique paralog list.
# And get the number of duplicate partners for each gene.

use strict;
use warnings;

open OUT, '>nfur/nothobranchius_furzeri_paralogs_with_counts.txt' or die $!;
open FH, 'nfur/nothobranchius_furzeri_paralogs.txt' or die $!;

#open OUT, '>aaul/aphyosemion_australe_paralogs_with_counts.txt' or die $!;
#open FH, 'aaul/aphyosemion_australe_paralogs.txt' or die $!;

#open OUT, '>olat/oryzias_latipes_paralogs_with_counts.txt' or die $!;
#open FH, 'olat/oryzias_latipes_paralogs.txt' or die $!;

#open OUT, '>drer/danio_rerio_paralogs_with_counts.txt' or die $!;
#open FH, 'drer/danio_rerio_paralogs.txt' or die $!;

#open OUT, '>alim/austrofundulus_limnaeus_paralogs_with_counts.txt' or die $!;
#open FH, 'alim/austrofundulus_limnaeus_paralogs.txt' or die $!;

# Read the parsed paralog file
my @file = <FH>;
my $head = shift @file;
chomp $head;
print OUT "$head\t";

my %Duplicates; # hash to store paralog pairs in each direction
my %Dupcounts; my %Dupnodes;

foreach (@file){ # Foreach paralog pair
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ((exists $Duplicates{$line[0]}{$line[1]}) || (exists $Duplicates{$line[1]}{$line[0]})){ # Check if it already exists in one or the other direction.
		
		die "pair $line[0] - $line[1] exists already\n";
	}
	else { # If not, increase the paralog counter by 1
		$Dupcounts{$line[0]} += 1;
		$Dupcounts{$line[1]} += 1;
		
		# And mark it checked
		push @{$Dupnodes{$line[0]}}, $line[2];
		push @{$Dupnodes{$line[1]}}, $line[2];
	}
	
}

#foreach (keys %Dupcounts){
	
#	print "$_\t$Dupcounts{$_}\t";
#	print join ',', @{$Dupnodes{$_}};
#	print "\n";
#}

# Print the unique paralogs in a new file with the partner count for each gene
print OUT "Total paralogs for Id1	Total paralogs for Id2	Sum total paralog pairs\n";
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

