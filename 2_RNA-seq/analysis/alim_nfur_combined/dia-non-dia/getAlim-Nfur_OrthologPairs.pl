# This is to get the corresponding paralog pairs between Nfur and Alim
use strict;
use warnings; # No warnings allowed

# Table to get proper gene names for Alim
#open FH, "/Volumes/Mybook_3/Other_organism_genomes/alim/Parsed_files/Annotation-table_Alim-ncbi100_withUTRs.txt" or die $!;
open ANNOT, "../Annotation-table_Alim-ncbi100_withUTRs.txt" or die $!;
my %Ids;

foreach (<ANNOT>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\"//g} @line;
	#print "$line[5]\t$line[0]\n";
	
	$Ids{$line[5]} = $line[0];
}

# Orthologs between alim and nfur
#open ORTH, "../Orthologues/Orthologues_nothobranchius_furzeri/nothobranchius_furzeri__v__austrofundulus_limnaeus.tsv" or die $!;
open ORTH, "../nothobranchius_furzeri__v__austrofundulus_limnaeus.tsv" or die $!;

# Ortholog files ----------------------------------------------------
my @orth = <ORTH>;
shift @orth;

my %Orthologs;

foreach (@orth){
	
	my @line = split "\t", $_;
	map {$_ =~s/\r\n|\"//g} @line;
	$line[1] =~s/NFU\|//g;
	$line[2] =~s/alimnaeus\|//g;
	
	if ($line[1] !~ /,/g && $line[2] !~ /,/g){ # Only one to one orthologs
		 #print "$line[1]\t$Ids{$line[2]}\n";
		 $Orthologs{$line[1]} = $Ids{$line[2]};
	}
}


# Paralog in other species --------------------------------------------
open PARA, "../../alim_dia_nondia/alimnaeus_DuplicatesClassification_All_le20_Reordered.txt" or die $!;
my @alim = <PARA>;
shift @alim;

my %AlimPara;

foreach (@alim){
	
	my @line = split "\t", $_;
	map {$_ =~s/\n|\"//g} @line;
 	$AlimPara{$line[0]}{$line[1]} = \@line;
}

open OUT, ">nfurzeri_alimnaeus_paralogs.txt" or die $!;
print OUT "nfur.gene1	nfur.gene2	nfur.duptime	nfur.category	";
print OUT "alim.gene1	alim.gene2	alim.duptime	alim.category	";
print OUT "median.Development.1	median.Diapause.1	median.FC.Diapause.1	median.FC.padj.1	median.Development.2	median.Diapause.2	median.FC.Diapause.2	median.FC.padj.2\t";
print OUT "Dia.log2FC.alim.1	Dia.padj.alim.1	median.Development.alim.1	median.Diapause.alim.1	Dia.log2FC.alim.2	Dia.padj.alim.2	median.Development.alim.2	median.Diapause.alim.2\n";

# Read Nfur paralogs and get the orthologs pairs in other species -----------------------
open NFUR, "../../nfur_diapause_CK/4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!;

my @nfur = <NFUR>;
shift @nfur;

foreach (@nfur){ # Foreach paralog pair in Alim
	
	my @line = split "\t", $_;
	map {$_ =~s/\n|\"//g} @line;
 	#$NfurPara{$line[0]}{$line[1]} = \@line;
 	
 	
 	if ((exists $Orthologs{$line[0]}) && (exists $Orthologs{$line[1]})){ # If both the gene in a pair have ortholog in Alim
 		
 		if (exists $AlimPara{$Orthologs{$line[0]}}{$Orthologs{$line[1]}}){	# And if they are a paralog in Alim in one direction
 			
 			my @alim = @{$AlimPara{$Orthologs{$line[0]}}{$Orthologs{$line[1]}}}; # Get Alim values
 			print OUT "$line[0]	$line[1]	$line[2]	$line[18]	$Orthologs{$line[0]}	$Orthologs{$line[1]}	$alim[2]	$alim[11]\t";
 			print OUT "$line[36]	$line[37]	$line[38]	$line[39]	$line[55]	$line[56]	$line[57]	$line[58]\t";
 			print OUT "$alim[12]	$alim[13]	$alim[22]	$alim[23]	$alim[24]	$alim[25]	$alim[34]	$alim[35]\n";
 			
 		}	
 		elsif (exists $AlimPara{$Orthologs{$line[1]}}{$Orthologs{$line[0]}}){ # If they are paralog in the other direction
 			my @alim = @{$AlimPara{$Orthologs{$line[1]}}{$Orthologs{$line[0]}}};
 			print OUT "$line[0]	$line[1]	$line[2]	$line[18]	$Orthologs{$line[1]}	$Orthologs{$line[0]}	$alim[2]	$alim[11]\t"; # Get ALim values
 			print OUT "$line[36]	$line[37]	$line[38]	$line[39]	$line[55]	$line[56]	$line[57]	$line[58]\t";
 			print OUT "$alim[24]	$alim[25]	$alim[34]	$alim[35]	$alim[12]	$alim[13]	$alim[22]	$alim[23]\n";
 		}
 	} 	
}

print `date`;



