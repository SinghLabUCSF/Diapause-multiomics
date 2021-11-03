# Read filtered diapause file with status and print the gene that is up first and down second for NeoF categories
use strict;
use warnings;

#######################################################################################################################

# ************* 71 vertebrates orthofider - main one for the paper *******************
open DUP, '4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!; # Outfile

# ************* 31 vertebrates orthofider ******************* ---------------------------------------------
#open DUP, '4_Paralog_classification/2_Orthofinder_31spp/nfurzeri_DuplicatesClassification_All_le20.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/2_Orthofinder_31spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!; # Outfile

# ************* 13 vertebrates orthofider ******************* ---------------------------------------------
#open DUP, '4_Paralog_classification/3_Orthofinder_13spp/nfurzeri_DuplicatesClassification_All_le20.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/3_Orthofinder_13spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!; # Outfile

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
		
	if ($line[18] eq 'NeoF'){ # If the category is neofunctionalized
		
		# If first one is up in diapause, it's all fine so just print
		if ((exists $UpInDiapauseCore{$line[0]}) && (exists $DownInDiapause{$line[1]})){

			print OUT join ("\t", @line), "\n";
			
		}
		# If second one is up in diapause - REVERSE THE ORFDER OF COLUMNS CORESPONDING TO GENE 1 AND GENE 2
		elsif ((exists $UpInDiapauseCore{$line[1]}) && (exists $DownInDiapause{$line[0]})){
			
			print OUT "$line[1]\t$line[0]\t";
			print OUT join ("\t", @line[2..7]), "\t";
			print OUT "$line[9]\t$line[8]\t";
			print OUT join ("\t", @line[10..20]), "\t";
			print OUT join ("\t", @line[40..58]), "\t";
			print OUT join ("\t", @line[21..39]), "\n";												
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
	elsif (($line[18] eq 'Both up') || ($line[18] eq 'Both down') || ($line[18] eq 'Unclassified')){ # Other categories, put the gene with higher FC in diapause first
		
		if ($line[38] > $line[57]){ # If the first gene FC is higher keep as is
				print OUT join ("\t", @line), "\n";		
		}
		else { # Else reverse the order
			print OUT "$line[1]\t$line[0]\t";
			print OUT join ("\t", @line[2..7]), "\t";
			print OUT "$line[9]\t$line[8]\t";
			print OUT join ("\t", @line[10..20]), "\t";
			print OUT join ("\t", @line[40..58]), "\t";
			print OUT join ("\t", @line[21..39]), "\n";														
		}
	}
	else {die "Category name *$line[18]* not found";} # Die if name does not matches
}

print `date`;
