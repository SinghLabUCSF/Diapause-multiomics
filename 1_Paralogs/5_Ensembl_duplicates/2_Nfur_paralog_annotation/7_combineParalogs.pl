# Combine paralogs from gene families and Ensmebl into a single final file
# Since I have checked earlier there is no duplicate node in the two files - I just filter and clean it
# I merge lineage specific duplciates in N. furzeri duplicates

use strict;
use warnings;

open OUT, ">8_final_nfur_duplicates/Final_nfurzeri_paralogs.txt" or die $!;
print OUT "Id1	Id2	Node	Source	Fish support (Ensembl node)\n";

# Ensembl duplication nodes -----------------------------------------------------------
open FH, "7_get_duplication_time_from_families/nfur_duplicates_from_families.txt" or die $!;
my @fam = <FH>;
shift @fam;

my %Families;

foreach (@fam){
		
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	print OUT "$line[0]\t$line[1]\t$line[5]\tFamilies\t0\n";
}


# Ensembl duplication nodes -----------------------------------------------------------
open ENS, "5_Combined_nodes/nfurzeri_consensus_dup_nodes.txt" or die $!;

foreach (<ENS>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	print OUT "$line[0]	$line[1]	";
	
	# Merge the following lineage specific nodes in N. furzeri
	if ($line[2] eq 'Otophysa' || 
		$line[2] eq 'Eupercaria' || 
		$line[2] eq 'Tetraodontidae' || 
		$line[2] eq 'Gasterosteus aculeatus' || 
		$line[2] eq 'Danio rerio' || 
		$line[2] eq 'Oryzias latipes' || 
		$line[2] eq 'Takifugu rubripes' || 
		$line[2] eq 'Tetraodon nigroviridis'){
		
		print OUT "Nothobranchius furzeri (lineage specific)	";
	}
	else {
		print OUT "$line[2]	";
	}
	print OUT "Ensembl\t$line[3]\n";
	
}

print `date`;
