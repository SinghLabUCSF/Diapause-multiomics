use strict;
use warnings;
use List::Util qw( min max );

# Outfiles from most relaxed to most stringent criteria. Numbers are from my analysis

# All: 46201 pairs

# All pairs after filtering more than 20 paralog partners 
open R0, ">List_All_20paralogs.txt"; # 20091 *************************************************** LIST 1
open R0A, ">List_All_20paralogs_best_matches.txt";

# >=40% similarity (33883 pairs)
open R1, ">List_Relaxed_1.txt";
# >=50% similarity (25137 pairs)
open R2, ">List_Relaxed_2.txt";
# >=50% identity <=25% gaps (11730 pairs)
open R3, ">List_Relaxed_3.txt";

# Highest similarity in family (10155 pairs)
open I1, ">List_Intermediate_1.txt"; # ********************************************************* LIST 2
# Highest similarity in family with >=50% similarity 
open I2, ">List_Intermediate_2.txt";

# Highest similarity in family non-redundant
open I3, ">List_Intermediate_3.txt";
# Highest similarity in family non-redundant >50% similarity
open I4, ">List_Intermediate_4.txt";

# Only duplicated once (2228 pairs) ************************************************************* LIST 3
open S1, ">List_Strict_1.txt";
# Only duplicated once with >50% similarity (1470 pairs)
open S2, ">List_Strict_2.txt";

open FH, "nfur_paralogs_KaKs_31spp_highestSimilarity.txt" or die $!;
my @file = <FH>;
my $head = shift @file;

print R1 "$head"; print R2 "$head"; print R3 "$head";
print I1 "$head"; print I2 "$head"; print I3 "$head"; print I4 "$head"; 
print S1 "$head"; print S2 "$head"; print R0 "$head"; print R0A "$head";

my %Checked;
my %Dupnodes; 
my %AllPairs;
my %LargeFamilies;

# Read all duplicate file and make 2 hashes: one for All pairs and one with only Ks in both directions
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	
	if ($line[10] <= 20){print R0 join ("\t", @line),"\n";} # All paralogs that have no more than 20 partners

	if ($line[10] <= 20){print R0A join ("\t", @line),"\n";} # All paralogs that have no more than 20 partners, with less than 20 in large families
	else {$LargeFamilies{$line[0]}{$line[1]}  = $line[12];
	      $LargeFamilies{$line[0]}{$line[1]}  = $line[12];
	}
		
	if ($line[10] == 2){print S1 join ("\t", @line),"\n";} # Pairs only duplicated once in the genome
	if ($line[10] == 2 && $line[12] >= 50){print S2 join ("\t", @line),"\n";} # Pairs only duplicated once in the genome with >50% similarity
		
	if ($line[12] >= 40){print R1 join ("\t", @line),"\n";} # Similarity >= 40%
	if ($line[12] >= 50){print R2 join ("\t", @line),"\n";} # Similarity >= 40%
	if ($line[11] >= 50 && 	$line[13] <= 25){print R3 join ("\t", @line),"\n";}	# Identity >= 50% & Gaps <= 25%

	$AllPairs{$line[0]}{$line[1]}  = \@line; # Similarity column
	$AllPairs{$line[1]}{$line[0]}  = \@line;
		
	# This is to get the highest similarity pairs for each family
	if ($line[17] == 1){ # 1 here means it is the best similarity pair among all of others
		
		print I1 join ("\t", @line),"\n"; # Keep only 1 highest similar pair for multi gene familes. Each gene can have at max 2 pairs
		if ($line[12] >= 50){print I2 join ("\t", @line),"\n";} # Same as above with >50% similarity
		
		# I want to further remove the ones that are redundant which I will do later
		#print "$line[0]\t$line[1]\n";
		$Dupnodes{$line[0]}{$line[1]}  = $line[12]; # Similarity column
		$Dupnodes{$line[1]}{$line[0]}  = $line[12];
		
	}
}

# Select minimum value pair by removing the redundant duplicates by taking maximum similarity among them
foreach my $id1 (keys %Dupnodes){
	
	if (scalar keys %{$Dupnodes{$id1}} > 1){
	
		my $highest = max values %{$Dupnodes{$id1}};
		#print join ("\t", $id1, scalar keys %{$Dupnodes{$id1}}, $highest), "\n";
		
		foreach my $id2 (keys %{$Dupnodes{$id1}}){
			
			if ($highest == $Dupnodes{$id1}{$id2}){
				
				if ((not exists $Checked{$id1}) && (not exists $Checked{$id2})){
					print I3 join ("\t", @{$AllPairs{$id1}{$id2}}), "\n";
					if (${$AllPairs{$id1}{$id2}}[12] >= 50){print I4 join ("\t", @{$AllPairs{$id1}{$id2}}), "\n";}
					$Checked{$id1} = '';
					$Checked{$id2} = '';
				}
			}
		}
	}
	else {
		foreach my $id2 (keys %{$Dupnodes{$id1}}){
			#print join ("\t", "\t->", $id1, scalar keys %{$Dupnodes{$id1}}), "\n";
			if ((not exists $Checked{$id1}) && (not exists $Checked{$id2})){
				print I3 join ("\t", @{$AllPairs{$id1}{$id2}}), "\n";
				if (${$AllPairs{$id1}{$id2}}[12] >= 50){print I4 join ("\t", @{$AllPairs{$id1}{$id2}}), "\n";}
				$Checked{$id1} = '';
				$Checked{$id2} = '';
			}
		}
	}
}

# For the families with more than 20 partners chose the most similar 20
my %Done;
foreach my $id1 (keys %LargeFamilies){

	my $i = 0;
	foreach my $id2 (sort {$LargeFamilies{$id1}{$b} cmp $LargeFamilies{$id1}{$a}} keys %{$LargeFamilies{$id1}}){

		#print "$id1\t$id2\t$LargeFamilies{$id1}{$id2}\n";
		print R0A join ("\t", @{$AllPairs{$id1}{$id2}}), "\n";
		$Done{$id1} += 1;
		$Done{$id2} += 1;

		$i++;
		if ($Done{$id1} == 10 || $Done{$id2} == 10){last;}
	}
}

print `date`;
