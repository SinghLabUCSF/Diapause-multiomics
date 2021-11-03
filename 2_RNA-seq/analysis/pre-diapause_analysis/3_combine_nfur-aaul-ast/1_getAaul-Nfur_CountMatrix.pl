# Get the RAW counts from Ast, Aaul file
# Get the RAW counts from Nfur file
# Make a single file that has Nfur symbols for orthologs, discard the others
# I will normalize it and use the normalized counts in later steps
use strict;
use warnings; # No warnings allowed

# Orthologs between alim and nfur
open ORTH, "nothobranchius_furzeri__v__aphyosemion_australe.tsv" or die $!;

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
		 #print "$line[1]\t$line[2]\n";
		 $Orthologs{$line[1]} = $line[2];
	}
}


# Counts in other species --------------------------------------------
open PARA, "../2_aaul_ast/Counts_development_aaul-ast_young_KV.csv" or die $!;
my @alim = <PARA>;
my $head = shift @alim;
chomp $head;
#print "$head";

my %AlimPara;

foreach (@alim){
	
	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
	#print "$line[0]\n";
 	$AlimPara{$line[0]} = [@line[1..11]];
}

open OUT, ">nfur_aaul-ast_ortholog_count_matrix.csv" or die $!;

# Read Nfur count matrix and add other columns if they exist
open NFUR, "../1_nfur/Counts_nfur_young.csv" or die $!;

my @nfur = <NFUR>;
my $head2 = shift @nfur;
chomp $head2;
print OUT "$head2\,$head\n";

foreach (@nfur){
	
	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
 	
 	#print  "$line[0]\n";
 	if (exists $Orthologs{$line[0]}){ # If it has ortholog 
		
		#print "$Orthologs{$line[0]}\n";
		if (exists $AlimPara{$Orthologs{$line[0]}}){ # and ortholog is expressed
			print OUT "\"$line[0]\",";
			print OUT (join ',', @line[1..4]),',';
			print OUT (join ',', @{$AlimPara{$Orthologs{$line[0]}}}),"\n";
		}
 	} 	
}

print `date`;



