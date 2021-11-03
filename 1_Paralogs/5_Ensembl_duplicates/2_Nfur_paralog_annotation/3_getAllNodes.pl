#
# Adapted from Ohnolog V2 Ohnologs/Synteny_All_2016_03_09/2_Paralogs/2_getAllNodes.pl
#
# This is to reconcile the duplication node from different ensembl versions
# I consider 7 (then) latest versions v80 to 86
# I decide final version from majority rule.
# Run it one by one for each organism and make sure there are no errors and warnings.
#
use strict;
use warnings;

local $| = 1;

my %Symbols;
my %Duplicates;
my %AllNodes;

my $organism = 'olatipes';

open OUT, ">5_Combined_nodes\/$organism\_CombinedDuplicates_Ens80-86.txt" or die $!;
open NODES, ">5_Combined_nodes\/AllNodes_Ens80-86_$organism\.txt" or die $!;
print OUT "Id1	Id2	Symbol1	Symbol2	Reconciled node	";

my $fileCount = 0;

foreach (<4_get_duplication_time_ensembl/nfurzeri-paralogs_$organism\_*.txt>) {
	
	print "$_\n";
	$_ =~/4_get_duplication_time_ensembl\/nfurzeri-paralogs_$organism\_(.+)\.txt/g;
	my $version = $1;
	print OUT "$version\t";
	
	open FH,"$_" or die $!;
	my @file = <FH>;
	shift @file;
	close ($_);
		
	foreach (@file){

		my @line = split "\t" or die $!;
		map {$_=~s/\n//g} @line;
		
		#$Symbols{$line[0]} = $line[1];

		if ($line[0] ne '' && $line[1] ne '' &&  $line[2] ne '' && $line[0] ne 'NA' && $line[1] ne 'NA' &&  $line[2] ne 'NA'){ # Its duplicate timing and id should not be null or NA
			
			if ((not exists $Duplicates{$line[0]."\t".$line[1]}) && (not exists $Duplicates{$line[1]."\t".$line[0]})){ # check if tha pair exists in any direction
				$Duplicates{$line[0]."\t".$line[1]} = ['', '', '', '', '', '', '']; # initialize the hash having an array of 7 for each ensembl version from 80 to 86
				${$Duplicates{$line[0]."\t".$line[1]}}[$fileCount] = $line[2];
			}
			else {
				# Here if a pair exists in both directions, the node will be replaced with the one that appear later. Although it shouldn't matter as
				# in the same version a pair cannot be associated with 2 nodes. Both must be same.
				if (exists $Duplicates{$line[0]."\t".$line[1]}){${$Duplicates{$line[0]."\t".$line[1]}}[$fileCount] = $line[2];} # push the value in the appropriate array
				if (exists $Duplicates{$line[1]."\t".$line[0]}){${$Duplicates{$line[1]."\t".$line[0]}}[$fileCount] = $line[2];} # push the value in the appropriate array
			}
			
			# push the node in nodes array
			$AllNodes{$line[2]} += 1;
		}
	}
	$fileCount++;		
}
print OUT "\n";

# For all duplication nodes - find out the final class here based on majority rule
foreach (keys %Duplicates){
	
	# Print ids
	print OUT "$_\t";
	my ($id1, $id2) = split "\t", $_;
	
	# Print symbols if exists
	if (exists $Symbols{$id1}){print OUT "$Symbols{$id1}\t";}
	else {print OUT "\t";}
	
	if (exists $Symbols{$id2}){print OUT "$Symbols{$id2}\t";}
	else {print OUT "\t";}
	
	my $count = 0;
	my $nonBlank = 0;
	
	# initialize variables to decide if the variable is old, intermediate or recent
	my %category;
	
	# Populate hash with counts
	foreach (@{$Duplicates{$_}}){
		
		if ($_ ne ''){			
			$category{$_} += 1;
		}
	}
	
	my @sortedvalues;
	my @sortedkeys;
	foreach (sort {$category{$b} cmp $category{$a}} keys %category){
		
		#print "$_\t$category{$_}\n";
		push @sortedvalues, $category{$_};
		push @sortedkeys, $_;
	}

	if (scalar @sortedvalues == 1){
		print OUT "$sortedkeys[0]";
	}
	elsif ((defined $sortedvalues[1]) && $sortedvalues[0] > $sortedvalues[1]){
		print OUT "$sortedkeys[0]";
		
	}
	elsif ((defined $sortedvalues[1]) && $sortedvalues[0] == $sortedvalues[1]){		
		print OUT "${$Duplicates{$_}}[5]"; # 
	}
	else {
		print OUT "";
	}

	print OUT "\t";
	# Print the duplicate timing for different versions
	print OUT (join "\t", @{$Duplicates{$_}});	
	print OUT "\n";
	
}


# print the nodes in the nodes file -- check this file properly and see if all the nodes are covered in this script
foreach (keys %AllNodes){
	
	print NODES "$_\t$AllNodes{$_}\n";
}

print `date`;
