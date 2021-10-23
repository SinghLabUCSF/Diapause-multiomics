# For lipids with multiple fatty acids, filter them based on mscore and tscore
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use List::Util qw(max);

# Open out file
open OUT, ">output/embryo_lipids_filtered_3.txt" or die $!;

# Read lipid file
open LIPID, "output/embryo_lipids_filtered_2.txt" or die $!;
my @lip = <LIPID>;
my $head = shift @lip;
print OUT "$head";

my %AllLipids; my %LipidCounts;

# For each lipid 
foreach (@lip){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
		# For rows with multiple FA, separate and see if they are non-unique. If they are the same print them. Else check how m-score and t-score are and take the best.
		# Because this is filtered file, all of them must be the same class
		
		print OUT join ("\t", @line[0..40]), "\t"; # Print everything until ion
		
		my @fattyacids = split '\|', $line[41];
		map {$_ =~s/\s+//g} @fattyacids;
		#print "@fattyacids\t";
				
		# Get m-score and t-score in arrays
		my @mscore = split '\|', $line[46];
		my @tscore = split '\|', $line[47];
		map {$_ =~s/\s+//g} @mscore;		
		map {$_ =~s/\s+//g} @tscore;
		
		# get max m-score and t-score
		my $maxM = max(@mscore);
		my $maxT = max(@tscore);
	
		if (scalar @fattyacids > 1){
			
#			print "@fattyacids\t@mscore\t@tscore\t*$maxM*\t*$maxT*\t"; # Must be a single class for all rows
			
			# check if there is more than 1 max m-score
			my $maxMct = grep {$_ == $maxM} @mscore;
#			print "$maxMct\t";
			
			if ($maxMct > 1){ # If 2FA have are more than 1 maximum M scores
				
				my $maxTct = grep {$_ == $maxT} @tscore; # Calculate t scores
#				print "$maxTct\t";
				
					if ($maxTct > 1){ # If 2FA have are more than 1 maximum T scores
						
						print "@fattyacids cannot be resolved - check!\n";
					}
					else {

						# get the index of max value and print that FA
						for (my $i = 0; $i <= $#tscore; $i++){
							if ($tscore[$i] eq $maxT){print OUT "$fattyacids[$i]\t"}
						} 												
					}
			}
			else {
				
				# get the index of max value and print that FA
				for (my $j = 0; $j <= $#mscore; $j++){
					if ($mscore[$j] eq $maxM){print OUT "$fattyacids[$j]\t"}
				} 												
			}
			
			my @class = split '\|', $line[42]; # print class
			print OUT "$class[0]\t";
			
			print OUT join ("\t", @line[43..47]), "\n"; # Print everything until ion			
						
#			print "\n";
		}
		else { # if there is just 1 FA,print as is
			print OUT join ("\t", @line[41..47]), "\n"; # Print everything until ion					
		}
}

print `date`;
