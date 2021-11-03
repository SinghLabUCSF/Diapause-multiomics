# Read filtered diapause file with status and print the gene that is up first and down second for NeoF categories
# *********** Different file for Ensembl because the column numbers are different
use strict;
use warnings;

#######################################################################################################################

# ************* Ensembl paralogs ******************* ---------------------------------------------
open DUP, '4_Paralog_classification/4_Ensembl/nfurzeri_DuplicatesClassification_All_le20.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
open OUT, ">4_Paralog_classification/4_Ensembl/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!; # Outfile

#######################################################################################################################
open DPUP, '2_Counts_DE/dpUp_relaxed.txt' or die $!; # Diapause up
open DPDOWN, '2_Counts_DE/dpDown_relaxed.txt' or die $!; # Diapause down

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
		
	if ($line[15] eq 'NeoF'){ # If the category is neofunctionalized
		
		# If first one is up in diapause, it's all fine so just print
		if ((exists $UpInDiapauseCore{$line[0]}) && (exists $DownInDiapause{$line[1]})){

			print OUT join ("\t", @line), "\n";
			
		}
		# If second one is up in diapause - REVERSE THE ORFDER OF COLUMNS CORESPONDING TO GENE 1 AND GENE 2
		elsif ((exists $UpInDiapauseCore{$line[1]}) && (exists $DownInDiapause{$line[0]})){
			
			print OUT "$line[1]\t$line[0]\t";
			print OUT join ("\t", @line[2..4]), "\t";
			print OUT "$line[6]\t$line[5]\t";
			print OUT join ("\t", @line[7..17]), "\t";
			print OUT join ("\t", @line[37..55]), "\t";
			print OUT join ("\t", @line[18..36]), "\n";												
		}
		# For other categories die as they should not be here
		elsif ((exists $UpInDiapauseCore{$line[0]}) &&\n (exists $UpInDiapauseCore{$line[1]})){		
			die "Both up at @line\n";
		}
		elsif ((exists $DownInDiapause{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			die "Both down at @line\n";
		}
		else {
			die "Unclassified at @line\n";
		}		
	}	
	elsif (($line[15] eq 'Both up') || ($line[15] eq 'Both down') || ($line[15] eq 'Unclassified')){ # Other categories, put the gene with higher FC in diapause first
		
		if ($line[35] > $line[54]){ # If the first gene FC is higher keep as is
				print OUT join ("\t", @line), "\n";		
		}
		else { # Else reverse the order
			print OUT "$line[1]\t$line[0]\t";
			print OUT join ("\t", @line[2..4]), "\t";
			print OUT "$line[6]\t$line[5]\t";
			print OUT join ("\t", @line[7..17]), "\t";
			print OUT join ("\t", @line[37..55]), "\t";
			print OUT join ("\t", @line[18..36]), "\n";												
		}
	}
	else {die "Category name *$line[15]* not found";} # Die if name does not matches
}

print `date`;
