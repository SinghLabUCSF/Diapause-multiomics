# Combine duplication node, WGD and diapause up genes into a single file.
# This one will read the genes up and down in diapause, fold change file and ohnolog files adn combine them into one.
# *********** NO WARNINGS ALLOWED **************
use strict;
use warnings;

#######################################################################################################################
# ************* 71 vertebrates orthofider - main one for the paper *******************
open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20.txt" or die $!; # Outfile  

# ************* 31 vertebrates orthofider ******************* ---------------------------------------------
#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug11_31spp/List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/2_Orthofinder_31spp/nfurzeri_DuplicatesClassification_All_le20.txt" or die $!; # Outfile

# ************* 13 vertebrates orthofider ******************* ---------------------------------------------
#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug12_13spp/List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/3_Orthofinder_13spp/nfurzeri_DuplicatesClassification_All_le20.txt" or die $!; # Outfile 

# ************* Ensembl paralogs *******************  -----------------------------------------------------
#open DUP, '3_Paralog_datasets/5_Ensembl_duplicates/List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/4_Ensembl/nfurzeri_DuplicatesClassification_All_le20.txt" or die $!; # Outfile

# ********** Other criteria *************
#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/List_All_20paralogs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_StrictLists.txt" or die $!; # Outfile       

#open DUP, '5_orthofinder_Results_Aug04_71spp/List_Strict_1.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">nfurzeri_DuplicatesClassification_Only1DupEvent.txt" or die $!; # Outfile

#open DUP, '5_orthofinder_Results_Aug04_71spp/List_Intermediate_1.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">nfurzeri_DuplicatesClassification_OnlyMostSimilar.txt" or die $!; # Outfile

# ************* Randomized paralogs *******************  -----------------------------------------------------
#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/List_All_20paralogs_ShuffledPairs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_ShuffledPairs.txt" or die $!; # Outfile

#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/nfur_paralogs_RandomizedPairs_3.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_RandomizedPairs_2.txt" or die $!; # Outfile

#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/nfur_paralogs_RandomizedPairs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_RandomizedPairs_StrictLists.txt" or die $!; # Outfile

#open DUP, '3_Paralog_datasets/5_orthofinder_Results_Aug04_71spp/nfur_paralogs_RandomizedSingletonPairs.txt' or die $!; # All filtered duplicates file with Ka, Ks etc
#open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_RandomizedSingletonPairs.txt" or die $!; # Outfile

#######################################################################################################################
open DPUP, '2_Counts_DE/dpUp_relaxed.txt' or die $!; # Diapause up. USED FOR PAPER
open DPDOWN, '2_Counts_DE/dpDown_relaxed.txt' or die $!; # Diapause down. USED FOR PAPER

#open DPUP, 'dpUp_intermediate.txt' or die $!; # Diapause up. for test
#open DPDOWN, 'dpDown_intermediate.txt' or die $!; # Diapause down. for test

#open DPUP, '2_Counts_DE/dpUp_strict.txt' or die $!; # Diapause up. for test
#open DPDOWN, '2_Counts_DE/dpDown_strict.txt' or die $!; # Diapause down. for test

open DPALL, '2_Counts_DE/Counts_DE_Combined_Final.csv' or die $!; # Diapause all

open OHNO, '../Datasets/nfurzeri_Criteria-[C]-Pairs_3R.txt' or die $!; # 3R Ohno
open OHNO2, '../Datasets/nfurzeri_Criteria-[C]-Pairs_2R.txt' or die $!; # 2R ohno
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
shift @dpall;
my %AllDE;

foreach (@dpall){
	
	my @line = split ',', $_;
	map {$_=~s/\n|\"//g} @line;
	#print "*$line[0]*\t$line[16]\t$line[17]\t$line[23]\t$line[24]\n";
	$AllDE{$line[0]} = [@line[15..29], $line[31], $line[32], $line[38], $line[39]]; # Remove NAs

}
print scalar keys %AllDE, "\n";
#print @{$AllDE{'LOC107377912'};
close (DPALL);



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


# All filtered duplicates file with Ka, Ks etc
my @dup = <DUP>;
my $head = shift @dup;
my @head = split "\t", $head;
map {$_=~s/\n//g} @head;

print OUT join ("\t", @head, "Category\t2R\t3R\tPreD1.1	PreD2.1	PreD3.1	D3d1.1	D3d2.1	D3d3.1	D6d1.1	D6d2.1	D6d3.1	D1m1.1	D1m2.1	D1m3.1	NonD1.1	NonD2.1	NonD3.1\tmedian.Development.1\tmedian.Diapause.1\tmedian.FC.Diapause.1\tmedian.FC.padj.1\tPreD1.2	PreD2.2	PreD3.2	D3d1.2	D3d2.2	D3d3.2	D6d1.2	D6d2.2	D6d3.2	D1m1.2	D1m2.2	D1m3.2	NonD1.2	NonD2.2	NonD3.2\tmedian.Development.2\tmedian.Diapause.2\tmedian.FC.Diapause.2\tmedian.FC.padj.2"), "\n";

my %Checked;

foreach (@dup){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
		
	if ((not exists $Checked{$line[0]}{$line[1]}) && (not exists $Checked{$line[1]}{$line[0]})){
		
		$Checked{$line[1]}{$line[0]} = '';
		$Checked{$line[0]}{$line[1]} = '';

		print OUT join ("\t", @line), "\t";	

		# One up in diapause and the other one in down
		if ((exists $UpInDiapauseCore{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			print OUT "NeoF\t";
		}
		# One up in diapause and the other one or down - reverse direction
		elsif ((exists $UpInDiapauseCore{$line[1]}) && (exists $DownInDiapause{$line[0]})){
			print OUT "NeoF\t"; # print in reverse so that 1st gene is the one that is up 
		}
		# Both up in diapause
		elsif ((exists $UpInDiapauseCore{$line[0]}) && (exists $UpInDiapauseCore{$line[1]})){
			print OUT "Both up\t";
		}
		# Both down in diapause
		elsif ((exists $DownInDiapause{$line[0]}) && (exists $DownInDiapause{$line[1]})){
			print OUT "Both down\t"; # print in reverse so that 1st gene is the one that is up 
		}
		# All other patterns
		else {
			print OUT "Unclassified\t";
			#print "$line[0]\t$line[1]\n"
		}

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
		
		if (exists $AllDE{$line[0]}){
			print OUT join ("\t", @{$AllDE{$line[0]}}),"\t";
		}
		else {print OUT "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";}
		
		if (exists $AllDE{$line[1]}){
			print OUT join ("\t", @{$AllDE{$line[1]}}),"\n";
		}
		else {print OUT "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";}
		
	}	
}

print `date`;
