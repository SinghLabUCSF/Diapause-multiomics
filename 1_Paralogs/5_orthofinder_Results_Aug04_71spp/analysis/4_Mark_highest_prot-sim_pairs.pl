# From the KaKs file, for the paralogs with more than one duplicate partners, mark the pairs with the lowest protein similarity. This will select a single best matching paralog pairs for each family.
# There will be warnings for NA values, which can be ignored.
use strict;
use warnings;
use List::Util qw( min max );

open OUT, ">nfur/nfur_paralogs_KaKs_71spp_highestSimilarity.txt" or die $!;

# Read the paralog file with Ka, Ks and protein similarity information.
open FH, "nfur/nfur_paralogs_KaKs_71spp.txt" or die $!;
my @file = <FH>;
my $head = shift @file;
my @head = split "\t", $head;
map {$_=~s/\n//g} @head;

my %Duplicates;
my %Dupnodes; 
my %AllPairs;
my %maxValues;

# Read all duplicate file and make 2 hashes: one for All pairs (unique) and one with only protein similarity in both directions (redundant)
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[2] ne ''){ # Must have a duplication node
		$AllPairs{$line[0]}{$line[1]}  = \@line;
	
		$Dupnodes{$line[0]}{$line[1]}  = $line[12]; # column 12 is protein Similarity
		$Dupnodes{$line[1]}{$line[0]}  = $line[12]; # column 12 is protein Similarity
	}
}

# From the similarity hash, select maximum similarity value pair and put it in another hash
foreach my $id1 (keys %Dupnodes){ # For each pair
	
	my $highest = max values %{$Dupnodes{$id1}}; # Select the maximum value

	foreach my $id2 (keys %{$Dupnodes{$id1}}){ 

		if ($highest == $Dupnodes{$id1}{$id2}){ # Get the pair with the max value
			
			if ((not exists $Duplicates{$id1}{$id2}) && (not exists $Duplicates{$id2}{$id1})){ # If it has not been checked already
				
				$Duplicates{$id1}{$id2} = ''; # mark it checked
				$Duplicates{$id2}{$id1} = '';
				
				$maxValues{$id1}{$id2} = ''; # Add it to a new hash that has maximum similarity value pairs
				$maxValues{$id2}{$id1} = '';
				
				# Uncomment this to print all pairs that are 
				#print "$id1	$id2	$Dupnodes{$id1}{$id2}\n";
			}
		}
	}
}

# Go through all pairs and put 1 as last colum if it's a minimum value multi pair, else 0
print OUT join ("\t", @head), "\tHighest similarity for multi-pairs\n";
foreach my $id1 (keys %AllPairs){
	
	foreach my $id2 (keys %{$AllPairs{$id1}}){
				
		print OUT join ("\t", @{$AllPairs{$id1}{$id2}}), "\t";
		if (exists $maxValues{$id1}{$id2}){	# For each pair if it exists in the max similarity pair, add 1 as a flag to select them later.
			print OUT "1\n";
		}
		else {
			print OUT "0\n";
		}
	}
}

print `date`;
