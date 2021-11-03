use strict;
use warnings; # No warnings allowed

# Table to get proper gene names for Alim
#open FH, "/Volumes/Mybook_3/Other_organism_genomes/alim/Parsed_files/Annotation-table_Alim-ncbi100_withUTRs.txt" or die $!;
open ANNOT, "Annotation-table_Alim-ncbi100_withUTRs.txt" or die $!;
my %Ids;

foreach (<ANNOT>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\"//g} @line;
	#print "$line[5]\t$line[0]\n";
	
	$Ids{$line[5]} = $line[0];
}

# Orthologs between alim and nfur
open ORTH, "nothobranchius_furzeri__v__austrofundulus_limnaeus.tsv" or die $!;

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
		 $Orthologs{$Ids{$line[2]}} = $line[1];
	}
}

open OUT, ">Normalized-Median-Counts_alim-longitudinal_nfurOrth.csv" or die $!;

# Counts in other species --------------------------------------------
open CTS, "Normalized-Median-Counts_alim-longitudinal.csv" or die $!;
my @alim = <CTS>;
my $head = shift @alim;
chomp $head;
print OUT "\"nfur.genes\",$head\n";

my %AlimPara;

foreach (@alim){
	
	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
	#print "$line[0]\n";
	
	if (exists $Orthologs{$line[0]}){
		
		print OUT "\"$Orthologs{$line[0]}\",\"$line[0]\",";
		print OUT join (',', @line[1..16]),"\n";
	}
}

print `date`;



