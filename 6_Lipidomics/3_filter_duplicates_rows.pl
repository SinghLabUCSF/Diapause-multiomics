# This will filter if there are multiple lipids in different rows by taking the one with highest intensity
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use List::Util qw(max);

# Open out file
open OUT, ">output/embryo_lipids_filtered_2.txt" or die $!;
open DISCARDED, ">output/embryo_lipids_discarded_duplicates_1.txt" or die $!;

# Read lipid file
open LIPID, "output/embryo_lipids_filtered_1.txt" or die $!;
my @lip = <LIPID>;
my $head = shift @lip;
chomp ($head);
print OUT "$head\n";
print DISCARDED "$head\tReason for filtering out\n";

my %AllLipids; my %LipidCounts;

# For each lipid 
foreach (@lip){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# For rows where multiple lipids are identified, separate the classes and take the class into a hash
	# Because this is filtered file, this must be the same classe

	my @classes = split '\|', $line[42]; # lipid class column
	map {$_ =~s/\s+//g} @classes;
	my @unique_classes = uniq @classes;
	#print "@unique_classes\n"; # Must be a single class for all rows
	
	# replace brackets from fatty acid column
	$line[41] =~s/\(|\)|\s+//g;
	#print "$line[41]\n";
	
	$LipidCounts{$unique_classes[0]}{$line[41]} += 1;
	push @{$AllLipids{$unique_classes[0]}{$line[41]}}, $_;
		
}

# Filter duplicates 
foreach my $class (keys %LipidCounts){ # foreach class
	
	foreach my $lipid (keys %{$LipidCounts{$class}}){ # foreach lipid
		
		#print "$class	$lipid	$LipidCounts{$class}{$lipid}\n";
		if ($LipidCounts{$class}{$lipid} == 1){ # if there is just one print it as is

			my @lipid_lines = @{$AllLipids{$class}{$lipid}};
			print OUT join "\n", @lipid_lines;
		}
		else { # If multiple lipids
			
			# Get intensities
			my @intensities;
			my @lipid_lines = @{$AllLipids{$class}{$lipid}};
			#print @lipid_lines, "\n";
			foreach (@lipid_lines){
				
				my @line = split "\t", $_;
				map {$_=~s/\n//g} @line;
				push @intensities, $line[4];
			}

			# get max intensity and print line corresponding to max intensity
			my $max = max(@intensities);
			
			foreach (@lipid_lines){
				
				my @line = split "\t", $_;
				map {$_=~s/\n//g} @line;

				if ($line[4] ==  $max){
					print OUT "$_";
				}
				else {
					print DISCARDED (join "\t", @line), "\tDuplicate lipid with low intensity\n";
				}
			}
		}
	}
}

print `date`;
