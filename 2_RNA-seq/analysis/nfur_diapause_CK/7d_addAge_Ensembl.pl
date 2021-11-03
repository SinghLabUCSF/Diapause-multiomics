# Foreach node add the age of paralog or it's duplciation time to divide them into 3 categories.
use strict;
use warnings;
use List::Util 'shuffle';

# Read the paralog file
open PARA, '4_Paralog_classification/4_Ensembl/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt' or die $!;
open OUT, ">4_Paralog_classification/4_Ensembl/nfurzeri_DuplicatesClassification_All_le20_Reordered_withAge.txt" or die $!;

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
'Euteleostomi' => 'Ancient vertebrate',
'Neopterygii' => 'Recent all fish',
'Clupeocephala' => 'Recent all fish',
'Teleosts' => 'Recent all fish',
'Otophysi' => 'Recent all fish',
'Acanthomorphata' => 'Recent all fish',
'Percomorphaceae' => 'Recent all fish',
'Ovalentaria' => 'Recent all fish',
'Atherinomorphae' => 'Recent all killifish',
'Aplocheiloidei' => 'Recent all killifish',
'Nothobranchiadae' => 'Recent all killifish',
'Nothobranchius furzeri' => 'Recent all killifish',
'Nothobranchius furzeri (lineage specific)' => 'Recent all killifish');

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
