use strict;
use warnings; # No warnings allowed

# Orthologs between alim and nfur
open ORTH, "nothobranchius_furzeri__v__oryzias_latipes.tsv" or die $!;

# Ortholog files ----------------------------------------------------
my @orth = <ORTH>;
shift @orth;

my %Orthologs;

foreach (@orth){
	
	my @line = split "\t", $_;
	map {$_ =~s/\r\n|\"//g} @line;
	$line[1] =~s/NFU\|//g;
	$line[2] =~s/(ENSORLG\d{11})\|//g;
	my $id = $1;
	
	if ($line[1] !~ /,/g && $line[2] !~ /,/g){ # Only one to one orthologs
		 print "$line[1]\t$id\n";
		 $Orthologs{$id} = $line[1];
	}
}

open OUT, ">Normalized-Median-Counts_medaka_all_nfurOrth.csv" or die $!;

# Counts in other species --------------------------------------------
open CTS, "Normalized-Median-Counts_medaka_all.csv" or die $!;
my @alim = <CTS>;
my $head = shift @alim;
chomp $head;
print OUT "$head\n";

my %AlimPara;

foreach (@alim){
	
	my @line = split "\,", $_;
	map {$_ =~s/\n|\"//g} @line;
	#print "$line[0]\n";
	
	if (exists $Orthologs{$line[0]}){
		
		print OUT "\"$Orthologs{$line[0]}\",";
		print OUT join (',', @line[1..6]),"\n";
	}
}

print `date`;

