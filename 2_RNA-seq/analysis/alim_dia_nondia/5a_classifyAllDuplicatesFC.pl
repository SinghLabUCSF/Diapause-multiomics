# Combine duplication node, WGD and diapause up genes
# *********** NO WARNINGS ALLOWED **************
use strict;
use warnings;

#######################################################################################################################

open DUP, 'List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc *************
open OUT, ">alimnaeus_DuplicatesClassification_All_le20.txt" or die $!; # Outfile             *************

#######################################################################################################################
open DPUP, 'dpUp.txt' or die $!; # Diapause up 
open DPDOWN, 'dpDown.txt' or die $!; # Diapause down

open DPALL, 'Counts_DE_Combined_Final_alim_dia-nondia.csv' or die $!; # Diapause all

#open OHNO, '/Volumes/Mybook_3/Ohnologs/Synteny_All_2016_03_09/7_FilterOhnologs/nfurzeri/nfurzeri_Criteria-[C]-Pairs_3R.txt' or die $!; # 3R Ohno
#open OHNO2, '/Volumes/Mybook_3/Ohnologs/Synteny_All_2016_03_09/7_FilterOhnologs/nfurzeri/nfurzeri_Criteria-[C]-Pairs_2R.txt' or die $!; # 2R ohno
#open FOXO, '/Volumes/Mybook_3/Motif_enrichment/results/foxo3MotifsAll.txt' or die $!; # File with foxo motifs # Not using for now

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

# Diapause DE all -------------------------------------------------------------
my @dpall = <DPALL>;
my $dupHead = shift @dpall;
my @dupHead = split "\,", $dupHead;
shift @dupHead;
map {$_=~s/\n//} @dupHead;

my %AllDE;

foreach (@dpall){
	
	my @line = split ',', $_;
	map {$_=~s/\n|\"//g} @line;
	#print "*$line[0]*\t$line[16]\t$line[17]\t$line[23]\t$line[24]\n";
	$AllDE{$line[0]} = [@line[1..12]]; # Remove NAs

}
print scalar keys %AllDE, "\n";
#print @{$AllDE{'LOC107377912'};
close (DPALL);


=cut
# 3R Ohno file ----------------------------------------------------------
my @ohno = <OHNO>;
shift @ohno;
my %Ohno3R;

foreach (@ohno){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Ohno3R{$line[0]}{$line[1]} = '';
	$Ohno3R{$line[1]}{$line[0]} = '';
}
close (OHNO);

# 2R Ohno file ------------------------------------------------------------
my @ohno2 = <OHNO2>;
shift @ohno2;
my %Ohno2R;

foreach (@ohno2){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Ohno2R{$line[0]}{$line[1]} = '';  
	$Ohno2R{$line[1]}{$line[0]} = '';
}
close (OHNO2);

# Has foxo motif
#my @fox = <FOXO>;
#shift @fox;
#my %Foxo;
#
#foreach (@fox){
#	
#	my @line = split "\t", $_;
#	map {$_=~s/\n//g} @line;
#	$Foxo{$line[0]} = '';
#}
=cut

# All filtered duplicates file with Ka, Ks etc

my @dup = <DUP>;
my $head = shift @dup;
my @head = split "\t", $head;
map {$_=~s/\n//g} @head;

print OUT join ("\t", @head, "Category"), "\t";
foreach (@dupHead){print OUT "$_\.1\t"}
foreach (@dupHead){print OUT "$_\.2\t"}
print OUT "\n";

my %Checked;

foreach (@dup){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
		
	if ((not exists $Checked{$line[0]}{$line[1]}) && (not exists $Checked{$line[1]}{$line[0]})){
		
		$Checked{$line[1]}{$line[0]} = '';
		$Checked{$line[0]}{$line[1]} = '';
		
		print OUT join ("\t", @line), "\t";
		
		# One up in diapause and the other one not changed or down
		if ((exists $UpInDiapauseCore{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			print OUT "NeoF\t";
		}
		# One up in diapause and the other one not changed or down - reverse direction
		elsif ((exists $UpInDiapauseCore{$line[1]}) && (exists $DownInDiapause{$line[0]})){
			print OUT "NeoF\t"; # print in reverse so that 1st gene is the one that is up 
		}
		# Both up in diapause
		elsif ((exists $UpInDiapauseCore{$line[0]}) && (exists $UpInDiapauseCore{$line[1]})){
			print OUT "Both up\t";
		}
		elsif ((exists $DownInDiapause{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			print OUT "Both down\t"; # print in reverse so that 1st gene is the one that is up 
		}
		else {
			print OUT "Unclassified\t";
			#print "$line[0]\t$line[1]\n"
		}
=cut
		# 2R Ohno or not
		if (exists $Ohno2R{$line[0]}{$line[1]}){print OUT "2R\t";}
		else {print OUT "\t";}
		# 3R Ohno or not afer excluding 2R ohno
		if ((not exists $Ohno2R{$line[0]}{$line[1]}) && (exists $Ohno3R{$line[0]}{$line[1]})){print OUT "3R\t";}
		else {print OUT "\t";}

#		# Foxo motif in promoter or not
#		if    ((exists $Foxo{$line[0]}) && (exists $Foxo{$line[1]})){print OUT "both\t";}
#		elsif ((exists $Foxo{$line[0]}) && (not exists $Foxo{$line[1]})){print OUT "Id1\t";}
#		elsif ((not exists $Foxo{$line[0]}) && (exists $Foxo{$line[1]})){print OUT "Id2\t";}
#		elsif ((not exists $Foxo{$line[0]}) && (not exists $Foxo{$line[1]})){print OUT "\t";}	
#		else {die "error with foxo";}		
=cut	
		if (exists $AllDE{$line[0]}){
			print OUT join ("\t", @{$AllDE{$line[0]}}),"\t";
		}
		else {print OUT "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";}
		
		if (exists $AllDE{$line[1]}){
			print OUT join ("\t", @{$AllDE{$line[1]}}),"\n";
		}
		else {print OUT "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";}
		
	}	
}

print `date`;
