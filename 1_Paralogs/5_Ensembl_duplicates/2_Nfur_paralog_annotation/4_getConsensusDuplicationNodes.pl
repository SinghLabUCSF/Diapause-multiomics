# Adapted from Ohnologs V2 
#
# Combine duplication time from all from all 5 fish 
# Get a file with all duplication timings -> by taking the consensus of all fish
#
# There are some warnings due to no node for some pairs (blank). these can be ignored.
# The nodes that have no consensus and are not there in the lastest version are blank. These will be ignored later.
# If there is no consensus on any of the node, I select one randomly. These are only 482 nodes out of ~59k, ~0.8% of all data so doesn't matter that much.
#
use strict;
use warnings;
use List::Util qw(max);
use List::UtilsBy qw(max_by);

# Initialize a fish array for printing
my @fish = ('drerio','gaculeatus','olatipes','tnigroviridis', 'trubripes');

# Outfile for Ohnolog runs
open OUT, ">5_Combined_nodes\/nfurzeri_all_dup_nodes.txt" or die $!;
print OUT "Id1	Id2	";
foreach (@fish){
	print OUT "$_\t";
}
print OUT "\n";

my %Paralogs;

foreach (<5_Combined_nodes\/*_CombinedDuplicates_Ens80-86.txt>){
	
	$_=~/.+\/(.+)_CombinedDuplicates_Ens80-86.txt/g;
	my $org = $1;
	#print "$org\n";	

#	if ($org eq 'drerio'){ # test for one organism

	print "$_\n";	
	open FH, "$_" or die $!;
	my @paralogs = <FH>;
	shift @paralogs;	
	
	foreach (@paralogs){
		
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;

		$Paralogs{$line[0]}{$line[1]}{$org} = $line[4]; # 4 is the reconciled node form different versions
		$Paralogs{$line[1]}{$line[0]}{$org} = $line[4];
	}

#	} # end test for 1 gene
	
}

# Print all nodes form the 5 fish in a file in both orientations
foreach my $key1 (keys %Paralogs){

	foreach my $key2 (keys %{$Paralogs{$key1}}){

		print OUT "$key1	$key2	";
		
		foreach (@fish){
		
			if (exists $Paralogs{$key1}{$key2}{$_}){print OUT "$Paralogs{$key1}{$key2}{$_}\t";}
			else {print OUT "\t";}
		}
		print OUT "\n";
	}
}
close(OUT);

# get consensus --------------------------------
open ALL, "5_Combined_nodes/nfurzeri_all_dup_nodes.txt" or die $!;
my @nodes = <ALL>;
shift @nodes;

open CON, ">5_Combined_nodes/nfurzeri_consensus_dup_nodes.txt" or die $!;
my %Duplicates;

foreach (@nodes){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ((not exists $Duplicates{$line[0]}{$line[1]}) && (not exists $Duplicates{$line[1]}{$line[0]})){
		
		$Duplicates{$line[0]}{$line[1]} = '';
		$Duplicates{$line[1]}{$line[0]} = '';
			
		print CON join ("\t", @line[0..1]),"\t";
	
		my %NodeCount;
		foreach (@line[2..6]){		
			$NodeCount{$_} += 1 if $_ ne '';
		}
		
		# If there are multiple maximums, one that comes first is selected, whcih is random and should average out.
		my $highest = max_by { $NodeCount{$_} } keys %NodeCount;
		
		print CON "$highest\t";
		print CON "$NodeCount{$highest}\n";
		
		# This is to print the node and animal support to see how many nodes have multiple maximums		
		print join ("\t", @line[0..1]),"\t";
		foreach (keys %NodeCount){
			print "$_\_$NodeCount{$_}\t";
		}
		print "\n";
	}
}
