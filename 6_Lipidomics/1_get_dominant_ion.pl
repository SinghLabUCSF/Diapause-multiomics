# Read the lipid file and make another one with most frequent ion for each lipid class
# I will filter the major ion file manually and use it in the next step. 
use strict;
use warnings;
use Statistics::Basic qw(:all);

# Files -----------------------
open FH, "data/embryo_lipids_all.txt";
open OUT, ">output/embryo_ion_counts.csv";
# -----------------------------

print OUT "lipid.class,ion,count,median,unique.lipids\n";

my @file = <FH>;
my $head = shift @file;

my %ClassIonCounts;
my %ClassIntensity;
my %UniqueLipids;

foreach (@file){ # Read all lipids file line by line
	
	my @line = split "\t", $_;
	if ($line[40] ne '') { # Column 40 in a 0 based index is lipid ion
		
		#print "$line[40]\n";
		
		my @multi = split '\|', $line[40]; # Split if there are multiple ions
		map {$_=~s/^\s+|\s+$//} @multi;
		#print "@multi\n";
		
		foreach (@multi){
			
			#print "$_\t";
			if ($_=~/(.+)\(.+\)([\-\+].+)/g){ # Separate lipid class and ion and count each one
				#print "$1\t$2\n";
				$ClassIonCounts{$1}{$2} += 1;
				push @{$ClassIntensity{$1}{$2}}, $line[4];

				if (scalar @multi == 1){$UniqueLipids{$1}{$2}{$multi[0]} = '';}
			}
		}		
	}
}

# Print
foreach my $class (keys %ClassIonCounts){

	foreach my $ion (keys %{$ClassIonCounts{$class}}){
		my $median = median(@{$ClassIntensity{$class}{$ion}});
		$median =~s/\,//g;
		print OUT "$class,\"$ion\",$ClassIonCounts{$class}{$ion},$median,";
		print OUT scalar keys %{$UniqueLipids{$class}{$ion}}, "\n";
	}
}

print `date`;

