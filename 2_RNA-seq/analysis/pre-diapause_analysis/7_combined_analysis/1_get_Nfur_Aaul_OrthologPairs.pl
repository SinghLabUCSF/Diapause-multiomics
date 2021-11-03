use strict;
use warnings; # No warnings allowed

# Normalized counts in nfur, Ast, Aaul KV
open NFURCT, "../3_combine_nfur-aaul-ast/Normalized-Median-Counts_nfur-aaul-ast.csv" or die $!;
my @nfurct = <NFURCT>;
shift @nfurct;
my %NfurAaulAst;

foreach (@nfurct){

	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
 	$NfurAaulAst{$line[0]} = "$line[1]\t$line[2]\t$line[3]"; # Young KV median
}

# Normalized counts in Alim
open ALIM, "../4_alim_dia_longitudinal/Normalized-Median-Counts_alim-longitudinal_nfurOrth.csv" or die $!;
my @alimnormct = <ALIM>;
my $alimcthead = shift @alimnormct;
my @alimcthead = split "\,", $alimcthead;
shift @alimcthead;
my %Alim;

foreach (@alimnormct){

	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
 	$Alim{$line[0]} = $line[9]; # 30C all somite stage median
}


# Normalized counts in Olat
open OLAT, "../5_medaka/Normalized-Median-Counts_medaka_all_nfurOrth.csv" or die $!;
my @olatct = <OLAT>;
shift @olatct;
my %Olat;

foreach (@olatct){

	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
 	$Olat{$line[0]} = $line[3]; # Stage16
}


# Normalized counts in Drer
open DRER, "../6_zebrafish/Normalized-Median-Counts_zebrafish_all_nfurOrth.csv" or die $!;
my @drct = <DRER>;
shift @drct;
my %Drer;

foreach (@drct){

	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
 	$Drer{$line[0]} = $line[3]; # 10-16hpf
}

open OUT, ">nfurzeri_ortho_paralogs_with_all_species.txt" or die $!;
print OUT "nfur.gene1	nfur.gene2	nfur.duptime	nfur.category	";
print OUT "nfur.kv.1\taaul.kv.1\tast.kv.1\t";
print OUT "nfur.kv.2\taaul.kv.2\tast.kv.2\t";
print OUT "alim.kv30.1\talim.kv30.2\t";
print OUT "olat.st16.1\tolat.st16.2\t";
print OUT "drer.10_16hpf.1\tdrer.10_16hpf.2";
print OUT "\n";

# Read Nfur paralogs and get the orthologs pairs in other species -----------------------
open NFUR, "../../nfur_diapause_CK/4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!;

my @nfur = <NFUR>;
shift @nfur;
my $necount = 0;
foreach (@nfur){
	
	my @line = split "\t", $_;
	map {$_ =~s/\n|\"//g} @line;
 	
	print OUT "$line[0]	$line[1]	$line[2]	$line[18]	";

	# Nfur Aaul Ast
	if ((exists $NfurAaulAst{$line[0]}) && exists ($NfurAaulAst{$line[1]})){
		print OUT "$NfurAaulAst{$line[0]}\t$NfurAaulAst{$line[1]}\t";
	}
	else {print OUT "\t\t\t\t\t\t";}
	
	# Alim
	if ((exists $Alim{$line[0]}) && exists ($Alim{$line[1]})){
		print OUT "$Alim{$line[0]}\t$Alim{$line[1]}\t";
	}	
	else {print OUT "\t\t";}
	
	# Olat
	if ((exists $Olat{$line[0]}) && exists ($Olat{$line[1]})){
		print OUT "$Olat{$line[0]}\t$Olat{$line[1]}\t";
	}
	else {print OUT "\t\t";}
	
	# Drer
	if ((exists $Drer{$line[0]}) && exists ($Drer{$line[1]})){
		print OUT "$Drer{$line[0]}\t$Drer{$line[1]}";
	}
	else {print OUT "\t";}
	
	print OUT "\n";
}

print "-> Nodes with different duplication time: $necount\n";
print `date`;



