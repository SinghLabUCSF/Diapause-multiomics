use strict;
use warnings;

# Read the paralog file
open PARA, 'alimnaeus_DuplicatesClassification_All_le20_Reordered.txt' or die $!;
open OUT, ">NeoF_Counts_All_le20.txt" or die $!;
print OUT "Ref_count	Ref_percent	NeoF_count	NeoF_percent	Ref_SUM	NeoF_SUM\n"; # PRINT header

my @file = <PARA>;
shift @file;

# All nodes for printing
my @all_nodes = ('Opisthokonta',
'Bilateria',
'Chordata',
'Vertebrata',
'Gnathostomata',
'Actinopterygii',
'Osteoglossocephalai',
'Clupeocephala',
'Euteleosteomorpha',
'Acanthomorphata_solderfish',
'Acanthomorphata_cod-sticklenback',
'Euacanthomorphacea',
'Percomorphaceae_roundgoby',
'Percomorphaceae_salarias',
'Ovalentaria',
'Atherinomorphae',
'Cyprinodontiformes',
'Aplocheiloidei',
'Austrofundulus_Kriptolebias',
'austrofundulus_limnaeus'
);

my %Vertebrate_nodes = ('Opisthokonta' => '',
'Bilateria' => '',
'Chordata' => '',
'Vertebrata' => '',
'Gnathostomata' => '');

my %Fish_nodes = ('Actinopterygii' => '',
'Osteoglossocephalai' => '',
'Clupeocephala' => '',
'Euteleosteomorpha' => '',
'Acanthomorphata_solderfish' => '',
'Acanthomorphata_cod-sticklenback' => '',
'Euacanthomorphacea' => '',
'Percomorphaceae_roundgoby' => '',
'Percomorphaceae_salarias' => '',
'Ovalentaria' => '',
'Atherinomorphae' => '');

my %Killifish_nodes = ('Cyprinodontiformes',
'Aplocheiloidei',
'Austrofundulus_Kriptolebias',
'austrofundulus_limnaeus'
);

my %AgeDiffCounts;
my %AgeCounts;

my %AllCounts;
my %DiffCounts;
#my $Count2R; my $Count2RDiff;
#my $Count3R; my $Count3RDiff;
#my $sum2r = 0; my $sum3r = 0;
	
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;

	# global criteria	
#	if ($line[7] == 2){ # With only highest similarity pair, more genes ~10k
#	if ($line[14] == 1 && $line[9] >= 50){ # With gene only duplicated once, less genes ~1200
	
		# All count
		$AllCounts{$line[2]} += 1;
	
		# differential count
		if ($line[11] eq 'NeoF'){
			$DiffCounts{$line[2]} += 1;
		}	

=cut
		# 2R count
		if ($line[19] eq '2R'){

			$Count2R += 1;
		
			if ($line[18] eq 'NeoF'){
				$Count2RDiff += 1;
			}
		}

		# 3R count
		if ($line[20] eq '3R'){
		
			$Count3R += 1;
		
			if ($line[18] eq 'NeoF'){
				$Count3RDiff += 1;
			}
		}
=cut		
		# Age counts
		if (exists $Vertebrate_nodes{$line[2]}){
			$AgeCounts{'Vertebrates'} += 1;
			if ($line[11] eq 'NeoF'){$AgeDiffCounts{'Vertebrates'} += 1;}
		}
		if (exists $Fish_nodes{$line[2]}){
			$AgeCounts{'Fish'} += 1;
			if ($line[11] eq 'NeoF'){$AgeDiffCounts{'Fish'} += 1;}
		}
		if (exists $Killifish_nodes{$line[2]}){
			$AgeCounts{'Killifish'} += 1;
			if ($line[11] eq 'NeoF'){$AgeDiffCounts{'Killifish'} += 1;}
		}
				
#	} # End global criteria
}

my $allsum = 0;
foreach (keys %AllCounts){$allsum = $allsum + $AllCounts{$_};}

my $diffsum = 0;
foreach (keys %DiffCounts){$diffsum = $diffsum + $DiffCounts{$_};}
	
# Print the counts
foreach (@all_nodes){
	
	print OUT "$_\t";
	if (exists $AllCounts{$_}){
		print OUT "$AllCounts{$_}\t";
		print OUT eval (($AllCounts{$_}*100)/$allsum), "\t";
	}
	else {print OUT "0\t0\t";}

	if (exists $DiffCounts{$_}){
		print OUT "$DiffCounts{$_}\t";
		print OUT eval (($DiffCounts{$_}*100)/$diffsum), "\t";
	}
	else {print OUT "0\t0\t";}
	
	print OUT "$allsum\t$diffsum\n";
		
}
#print OUT "\t$allsum\t\t$diffsum\n\n";
#print OUT "\n";
print "\tAll paralogs: $allsum\t\tNeoF paralog counts: $diffsum\n\n";


foreach (keys %AgeCounts){
	
	if (exists $AgeCounts{$_}){
		print OUT "$_\t$AgeCounts{$_}\t";
		print OUT eval (($AgeCounts{$_}*100)/$allsum), "\t";
	}
	else {print OUT "0\t0\n";}
	
	if (exists $AgeDiffCounts{$_}){
		print OUT "$AgeDiffCounts{$_}\t";
		print OUT eval (($AgeDiffCounts{$_}*100)/$diffsum), "\t";
	}
	else {print OUT "0\t0\t";}
	
	print OUT "$allsum\t$diffsum\n";
}
#print OUT "\n";

=cut
print OUT "2R WGD	$Count2R\t";
print OUT eval (($Count2R*100)/$allsum), "\t";
print OUT "$Count2RDiff\t";
print OUT eval (($Count2RDiff*100)/$diffsum), "\n";


print OUT "3R WGD	$Count3R\t";
print OUT eval (($Count3R*100)/$allsum), "\t";
print OUT "$Count3RDiff\t";
print OUT eval (($Count3RDiff*100)/$diffsum), "\n";
=cut

print `date`;