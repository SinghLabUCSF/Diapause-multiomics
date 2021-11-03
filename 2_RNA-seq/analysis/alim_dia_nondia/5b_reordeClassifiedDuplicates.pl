# Read filtered diapause file with status and print the gene that is up first and down second for NeoF categories
use strict;
use warnings;

#######################################################################################################################
open DUP, 'alimnaeus_DuplicatesClassification_All_le20.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
open OUT, ">alimnaeus_DuplicatesClassification_All_le20_Reordered.txt" or die $!; # Outfile

#######################################################################################################################
open DPUP, 'dpUp.txt' or die $!; # Diapause up
open DPDOWN, 'dpDown.txt' or die $!; # Diapause down

# Diapause up core -------------------------------------------------------------
my %UpInDiapauseCore;

foreach (<DPUP>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\"//g} @line;
	#print "$line[0]\n";
	$UpInDiapauseCore{$line[0]} = '' if ($line[0] ne 'NA'); # Remove NAs
}
close (DPUP);
print scalar keys %UpInDiapauseCore, "\n";

# Diapause down -------------------------------------------------------------
my %DownInDiapause;

foreach (<DPDOWN>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\"//g} @line;
	
	$DownInDiapause{$line[0]} = '' if ($line[0] ne 'NA'); # Remove NAs
}
close (DPDOWN);
print scalar keys %DownInDiapause, "\n";

# All filtered duplicates file with Ka, Ks etc

my @dup = <DUP>;
my $head = shift @dup;
my @head = split "\t", $head;
map {$_=~s/\n//g} @head;

print OUT join ("\t", @head), "\n";

my %Checked;

foreach (@dup){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
		
	if ($line[11] eq 'NeoF'){ # If the category is neofunctionalized
		
		# If first one is up in diapause, it's all fine so just print
		if ((exists $UpInDiapauseCore{$line[0]}) && (exists $DownInDiapause{$line[1]})){

			print OUT join ("\t", @line), "\n";
			
		}
		# If second one is up in diapause - REVERSE THE ORFDER OF COLUMNS CORESPONDING TO GENE 1 AND GENE 2
		elsif ((exists $UpInDiapauseCore{$line[1]}) && (exists $DownInDiapause{$line[0]})){
			
			print OUT "$line[1]\t$line[0]\t";
			print OUT join ("\t", @line[2..7]), "\t";
			print OUT "$line[9]\t$line[8]\t$line[10]\t$line[11]\t";
			print OUT join ("\t", @line[24..35]), "\t";
			print OUT join ("\t", @line[12..23]), "\n";
		}
		# For other categories die as they should not be here
		elsif ((exists $UpInDiapauseCore{$line[0]}) && (exists $UpInDiapauseCore{$line[1]})){
			die "Both up at @line\n";
		}
		elsif ((exists $DownInDiapause{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			die "Both down at @line\n";
		}
		else {
			die "Unclassified at @line\n";
		}		
	}	
	else {print OUT join ("\t", @line), "\n";} # For other categories just print as is
}

print `date`;
