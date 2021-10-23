# Use the major ion file and filter out the lipids with:
# No internal standard
# No major ion
# Belonging to multiple lipid classes
# Unannotated lipid
use strict;
use warnings;
use List::MoreUtils qw(uniq);

# Read major ion file
open IONS, "data/embryo_major_ions.txt" or die $!;
my @ions = <IONS>;
shift @ions;

my %MajorIons;

foreach (@ions){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$MajorIons{$line[0]}{$line[1]} = '';
}
#foreach (keys %MajorIons){print "$_\t$MajorIons{$_}\n";}

# Read normalization concentrations and make a hash for each tissue
open CONC, "data/standards_vol_for_normalization_Kevin.txt" or die $!;
my @conc = <CONC>;
shift @conc;

my %finalISconc;

foreach (@conc){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$finalISconc{"embryo"}{$line[1]} = $line[3];
	$finalISconc{"yolk"}{$line[1]} = $line[4];
	$finalISconc{"ovary"}{$line[1]} = $line[5];
}

# Open out file
open OUT, ">output/embryo_lipids_filtered_1.txt" or die $!;
open DISCARDED, ">output/embryo_lipids_discarded.txt" or die $!;

# Read lipid file
open LIPID, "data/embryo_lipids_all.txt" or die $!;
my @lip = <LIPID>;
my $head = shift @lip;
chomp ($head);
print OUT "$head\n";
print DISCARDED "$head\tReason for filtering out\n";

# For each lipid 
foreach (@lip){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# This will exclude all IS and unknown lipids
	if ($line[41] ne ''){
		
		# For rows where multiple lipids are identified, separate the classes and find how many unique classes are there
		#print "$line[42]\n";
		my @classes = split '\|', $line[42];
		map {$_ =~s/\s+//g} @classes;
		my @unique_classes = uniq @classes;
		#print "@unique_classes\n";


		# If only one unique class, then include them. 
		if (scalar @unique_classes == 1){
			#print "@unique_classes\n";
			
			# If the lipid class has IS
			if (exists $finalISconc{"embryo"}{$unique_classes[0]}){

				# Check if any of the ion is a major ion. If it is keep it, else discard. If theer are more than one lipid ions, it checks if any of them matches a major ion. I may need to filter them more later.
				my @lipidIons = split '\|', $line[40]; # lipid ions column
				map {$_ =~s/\s+//g} @lipidIons;

				my @ions; my $majorIonMatched = 0;

				foreach (@lipidIons){ # foreach ion
			
					$_=~/(.+)\(.+\)(.+)/g;
			
					if (exists $MajorIons{$1}{$2}){
						#print "$1\t$2\t";
						$majorIonMatched = 1;
					}
					#else {print " --> $1\t$2\t";}
				}
				#print "\t$majorIonMatched\t\n"
				
				if ($majorIonMatched == 1){
					print OUT join ("\t", @line),"\n";
				}
				else {print DISCARDED join ("\t", @line),"\tNot the major ion\n";}

			}
			else {print DISCARDED join ("\t", @line),"\tNo internal standard\n";} # Else filter out}
		}
		else {print DISCARDED join ("\t", @line),"\tMultiple lipid classes\n";} # Else filter out - ASK KEVIN what to do with these
			
	}
	else {print DISCARDED join ("\t", @line),"\tUnannotated lipid\n";} # Else filter out
}
print `date`;






