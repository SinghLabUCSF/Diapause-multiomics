# To generate a table of counts of paralogs belonging to each node and to the 3 age-categories. This table wil be used for stats in the paper.
use strict;
use warnings;

my $category = 'NeoF'; # NeoF, Unclassified
my $infile = '4_Paralog_classification/3_Orthofinder_71spp'; 
my $outfile = '5_Node_counts/3_Orthofinder_71spp';
my $out_suffix = 'le20'; # Only1DupEvent or le20 or OnlyMostSimilar. **** ALSO CHANGE THE PARALOG FILE BELOW.

# Read the paralog file
open PARA, "$infile/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!;   # Default 20 genes per family
#open PARA, "$infile/nfurzeri_DuplicatesClassification_OnlyMostSimilar.txt" or die $!;       # For a single most similar pair per family
#open PARA, "$infile/nfurzeri_DuplicatesClassification_Only1DupEvent.txt" or die $!;        # For genes duplicated only once


#open PARA, '4_Paralog_classification/1_Orthofinder_71spp/nfurzeri_DuplicatesClassification_All_le20_ShuffledPairs.txt' or die $!;

open OUT, ">$outfile/All_Counts_$category\_$out_suffix\.txt" or die $!;
print OUT "Ref_count	Ref_percent	$category\_count	$category\_percent\tRef\_SUM\t$category\_SUM\n"; # PRINT header

my @file = <PARA>;
shift @file;

# All nodes for printing
my @all_nodes = ('Opisthokonta',
'Bilateria',
'Chordata and above',
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
'Aplocheiloidei_pachypanchax',
'Nothobranchiidae_callopanchax',
'Nothobranchiidae_australe',
'Nothobranchius',
'Nothobranchius_furzeri');

my %Vertebrate_nodes = ('Opisthokonta' => '',
'Bilateria' => '',
'Chordata' => '',
'Chordata and above',
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

my %Killifish_nodes = ('Cyprinodontiformes' => '',
'Aplocheiloidei' => '',
'Aplocheiloidei_pachypanchax' => '',
'Nothobranchiidae_callopanchax' => '',
'Nothobranchiidae_australe' => '',
'Nothobranchius' => '',
'Nothobranchius_furzeri' => '');

my %AgeDiffCounts;
my %AgeCounts;

my %AllCounts;
my %DiffCounts;
my $Count2R = 0; my $Count2RDiff = 0;
my $Count3R = 0; my $Count3RDiff = 0;
my $sum2r = 0; my $sum3r = 0;
	
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;

	# global criteria	
#	if ($line[7] == 2){ # With only highest similarity pair, more genes ~10k
#	if ($line[14] == 1 && $line[9] >= 50){ # With gene only duplicated once, less genes ~1200
	
		# All count
		$AllCounts{$line[2]} += 1;
	
		# differential count
		if ($line[18] eq $category){
			$DiffCounts{$line[2]} += 1;
		}	

		# 2R count
		if ($line[19] eq '2R'){

			$Count2R += 1;
		
			if ($line[18] eq $category){
				$Count2RDiff += 1;
			}
		}

		# 3R count
		if ($line[20] eq '3R'){
		
			$Count3R += 1;
		
			if ($line[18] eq $category){
				$Count3RDiff += 1;
			}
		}
		
		# Age counts
		if (exists $Vertebrate_nodes{$line[2]}){
			$AgeCounts{'Vertebrates'} += 1;
			if ($line[18] eq $category){$AgeDiffCounts{'Vertebrates'} += 1;}
		}
		if (exists $Fish_nodes{$line[2]}){
			$AgeCounts{'Fish'} += 1;
			if ($line[18] eq $category){$AgeDiffCounts{'Fish'} += 1;}
		}
		if (exists $Killifish_nodes{$line[2]}){
			$AgeCounts{'Killifish'} += 1;
			if ($line[18] eq $category){$AgeDiffCounts{'Killifish'} += 1;}
		}
				
#	} # End global criteria
}

my $allsum = 0;
foreach (keys %AllCounts){$allsum = $allsum + $AllCounts{$_};}

my $diffsum = 0;
foreach (keys %DiffCounts){$diffsum = $diffsum + $DiffCounts{$_};}
	

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
print "\tAll paralogs: $allsum\t\t$category paralog counts: $diffsum\n\n";


foreach (keys %AgeCounts){
	
	print OUT "$_\t$AgeCounts{$_}\t";
	print OUT eval (($AgeCounts{$_}*100)/$allsum), "\t";
	
	print OUT "$AgeDiffCounts{$_}\t";
	print OUT eval (($AgeDiffCounts{$_}*100)/$diffsum), "\t";
	
	print OUT "$allsum\t$diffsum\n";
	
}
#print OUT "\n";

print OUT "2R_WGD	$Count2R\t";
print OUT eval (($Count2R*100)/$allsum), "\t";
print OUT "$Count2RDiff\t";
print OUT eval (($Count2RDiff*100)/$diffsum), "\t";

print OUT "$allsum\t$diffsum\n";

print OUT "3R_WGD	$Count3R\t";
print OUT eval (($Count3R*100)/$allsum), "\t";
print OUT "$Count3RDiff\t";
print OUT eval (($Count3RDiff*100)/$diffsum), "\t";

print OUT "$allsum\t$diffsum\n";

print `date`;
