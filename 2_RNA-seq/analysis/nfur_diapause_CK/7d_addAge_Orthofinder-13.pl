# Foreach node add the age of paralog or it's duplciation time to divide them into 3 categories.
use strict;
use warnings;
use List::Util 'shuffle';

# Read the paralog file
open PARA, '4_Paralog_classification/3_Orthofinder_13spp/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt' or die $!;
open OUT, ">4_Paralog_classification/3_Orthofinder_13spp/nfurzeri_DuplicatesClassification_All_le20_Reordered_withAge.txt" or die $!;

# Open the paralog file
my @file = <PARA>;
my $head = shift @file;
chomp $head;
print OUT $head, "\tAge category\n";

# All nodes for printing
my %Age = ('Chordata and above' => 'Ancient vertebrate',
'Vertebrata' => 'Ancient vertebrate',
'Gnathostomata,Euteleostomi,' => 'Ancient vertebrate',
'Actinopterygii,Neopterygii' => 'Recent all fish',
'Osteoglossocephalai,Clupeocephala' => 'Recent all fish',
'Euteleosteomorpha,Acanthomorphata' => 'Recent all fish',
'Euacanthomorphacea,Percomorphaceae,Ovalentaria' => 'Recent all fish',
'Atherinomorphae,Cyprinodontiformes' => 'Recent all killifish',
'Cyprinodontoidei' => 'Recent all killifish',
'Aplocheiloidei' => 'Recent all killifish',
'Nothobranchiidae' => 'Recent all killifish',
'nothobranchius_furzeri' => 'Recent all killifish');

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
