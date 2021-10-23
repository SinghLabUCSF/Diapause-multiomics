# Normalize lipids with respect to internal standards
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Statistics::Basic qw(:all);

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

# Read standards file and make a hash
open IS, "data/embryo_standards_filtered.txt" or die $!;
my @is = <IS>;
shift @is;

my %IS;

foreach (@is){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$IS{$line[42]} = [@line];
}

# Open out file
open OUT, ">output/embryo_lipids_normalized_IS.txt" or die $!;

# Read lipid file
open LIPID, "output/embryo_lipids_filtered_4.txt" or die $!;
my @lip = <LIPID>;
my $head = shift @lip;
print OUT $head;
my @head = split "\t", $head;

# For each lipid -- remember columns have samples
foreach (@lip){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# Because I have filetred the class there is just one here 
	#print "$line[42]\n";

	my $class = $line[42];
		
	print OUT join ("\t", @line[0..4]), "\t";
					
	for (my $i = 5; $i <= 39; $i++){ # For each sample
	
		my $conc = (($line[$i] * $finalISconc{"embryo"}{$class}) / ${$IS{$class}}[$i]);
		print OUT "$conc\t";
		
	}

	print OUT join ("\t", @line[40..47]), "\n";	# Last columns	
}
close (OUT);

print `date`;
