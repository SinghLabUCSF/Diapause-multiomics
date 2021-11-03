# Foreach node add the age of paralog or it's duplciation time to divide them into 3 categories.
use strict;
use warnings;
use List::Util 'shuffle';

# Read the paralog file
open PARA, '4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt' or die $!;
open OUT, ">4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_Reordered_withAge.txt" or die $!;

# Open the paralog file
my @file = <PARA>;
my $head = shift @file;
chomp $head;
print OUT $head, "\tAge category\n";

# All nodes for printing
my %Age = ('Opisthokonta' => 'Ancient vertebrate',
'Bilateria' => 'Ancient vertebrate',
'Chordata' => 'Ancient vertebrate',
'Vertebrata' => 'Ancient vertebrate',
'Gnathostomata' => 'Ancient vertebrate',
'Actinopterygii' => 'Recent all fish',
'Osteoglossocephalai' => 'Recent all fish',
'Clupeocephala' => 'Recent all fish',
'Euteleosteomorpha' => 'Recent all fish',
'Acanthomorphata_solderfish' => 'Recent all fish',
'Acanthomorphata_cod-sticklenback' => 'Recent all fish',
'Euacanthomorphacea' => 'Recent all fish',
'Percomorphaceae_roundgoby' => 'Recent all fish',
'Percomorphaceae_salarias' => 'Recent all fish',
'Ovalentaria' => 'Recent all fish',
'Atherinomorphae' => 'Recent all fish',
'Cyprinodontiformes' => 'Recent all killifish',
'Aplocheiloidei' => 'Recent all killifish',
'Aplocheiloidei_pachypanchax' => 'Recent all killifish',
'Nothobranchiidae_callopanchax' => 'Recent all killifish',
'Nothobranchiidae_australe' => 'Recent all killifish',
'Nothobranchius' => 'Recent all killifish',
'Nothobranchius_furzeri' => 'Recent all killifish');

# Foreach line
foreach (@file){
		
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	print OUT join ("\t", @line), "\t";
	
	if (exists $Age{$line[2]}){
		print OUT "$Age{$line[2]}\n";
	}
	else {die "$line[2] node not in Age hash!"}
}

print `date`;
