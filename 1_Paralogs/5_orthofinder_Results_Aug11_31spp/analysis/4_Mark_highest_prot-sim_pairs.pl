# From the KaKs file, the ones that have more than one duplicate partners, get the pairs that mark the pairs with the lowest Ka, Ks values
# There will be warnings for NA values, which can be ignored
use strict;
use warnings;
use List::Util qw( min max );

open OUT, ">nfur_paralogs_KaKs_31spp_highestSimilarity.txt" or die $!;

open FH, "nfur_paralogs_KaKs_31spp.txt" or die $!;
my @file = <FH>;
my $head = shift @file;
my @head = split "\t", $head;
map {$_=~s/\n//g} @head;

my %Duplicates;
my %Dupnodes; 
my %AllPairs;
my %maxValues;

# Read all duplicate file and make 2 hashes: one for All pairs and one with only Ks in both directions
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[2] ne ''){ # Ignore the ones with no node info because no consensus could be reached
		$AllPairs{$line[0]}{$line[1]}  = \@line;
	
		$Dupnodes{$line[0]}{$line[1]}  = $line[12]; # 12 is Similarity
		$Dupnodes{$line[1]}{$line[0]}  = $line[12]; # 12 is Similarity
	}
}

# From the similarity hash, select maximum similarity value pair and put it in another hash
foreach my $id1 (keys %Dupnodes){
	
	my $highest = max values %{$Dupnodes{$id1}};

	foreach my $id2 (keys %{$Dupnodes{$id1}}){

		if ($highest == $Dupnodes{$id1}{$id2}){
			
			if ((not exists $Duplicates{$id1}{$id2}) && (not exists $Duplicates{$id2}{$id1})){
				
				$Duplicates{$id1}{$id2} = '';
				$Duplicates{$id2}{$id1} = '';
				
				$maxValues{$id1}{$id2} = '';
				$maxValues{$id2}{$id1} = '';
				
				# Uncomment thsi to print all pairs that are 
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
		if (exists $maxValues{$id1}{$id2}){	
			print OUT "1\n";
		}
		else {
			print OUT "0\n";
		}
	}
}

print `date`;
