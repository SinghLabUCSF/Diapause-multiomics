use strict;
use warnings;
use List::MoreUtils qw(first_index);

# Variables 
my $organism = 'nothobranchius_furzeri'; # The name of organism as it appears in the duplicates.tsv file
my $pattern = 'NFU\|'; # pattern for gene id as it appears in the Orthogroups.tsv file to be replaced
my $outfile = 'nothobranchius_furzeri_singletons.txt'; # Outfile with all singletons

# Open outfiles
open OUT, ">$outfile" or die $!;
print OUT "gene\torthogroup\n";

# Read orthogroup count files
open CTS, "../Orthogroups/Orthogroups.GeneCount.tsv" or die $!; # Gene count to get singletons in the current species
my @counts = <CTS>;
my @head = split "\t", shift @counts;
#print "@head";
my $ctindex = first_index { $_ eq $organism } @head; # get the 0 based index of the organism column
#print "$ctindex\n";
#print "$counts[0]\n";
my %singletonIds;

foreach (@counts){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	#print "$line[0]\n";
	
	if ($line[$ctindex] == 1){ # If there is only a single gene for this orthogroup in selected organism
		
		#print "$line[0]\t$line[$ctindex]\n"
		$singletonIds{$line[0]} = @line;
	}		
}

# Read orthogroup Id files
open IDS, "../Orthogroups/Orthogroups.tsv" or die $!; # Gene count to get singletons in the current species
my @ids = <IDS>;
my @idhead = split "\t", shift @ids;
#print "@idhead\n";
my $idindex = first_index { $_ eq $organism } @idhead; # get the 0 based index of the organism column
#print "$idindex\n"; 


foreach (@ids){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	#print "$line[0]\n";
	
	if (exists $singletonIds{$line[0]}){ # If the orthogroup is the one with a single gene
		
		$line[$idindex] =~s/$pattern//g;
		print OUT "$line[0]\t$line[$idindex]\n"
		
	}		
}

=cut
# Read unassigned gene file
open USG, "../Orthogroups/Orthogroups_UnassignedGenes.tsv" or die $!; # Unassigned genes
my @unassigned = <USG>;
my @ushead = split "\t", shift @unassigned;
#print "@idhead\n";
my $usindex = first_index { $_ eq $organism } @ushead; # get the 0 based index of the organism column
#print "$idindex\n"; 


foreach (@unassigned){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	#print "$line[0]\n";
	if ($line[$usindex] ne ''){ # If theere are unassigned gene for the orgaanism
		print "$line[0]\t$line[$usindex]\n"
	}
}
=cut

print `date`;
