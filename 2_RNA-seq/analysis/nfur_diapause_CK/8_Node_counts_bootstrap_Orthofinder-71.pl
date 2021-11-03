# To do 10000 bootstrp with 50% of the data and count ancient, very ancient and recent specialized paralogs i each bootstrap.
# It will generate one file for each age category. One for each WGD and one that has everything.
# The file for each age categories are used for plotting.
use strict;
use warnings;
use List::Util 'shuffle';

# Input variables -----------------------------------------------------------------------------
my $category = 'NeoF'; # Unclassified or NeoF
my $infile = '4_Paralog_classification/1_Orthofinder_71spp'; 
my $outfile = '5_Node_counts/1_Orthofinder_71spp';
my $out_suffix = 'le20'; # Only1DupEvent or le20 or OnlyMostSimilar

# Read the paralog file ----------------------------------------------------------------------
open PARA, "$infile/nfurzeri_DuplicatesClassification_All_le20_Reordered.txt" or die $!;    # Default 20 genes per family
#open PARA, "$infile/nfurzeri_DuplicatesClassification_Only1DupEvent.txt" or die $!;        # For genes duplicated only once
#open PARA, "$infile/nfurzeri_DuplicatesClassification_OnlyMostSimilar.txt" or die $!;      # For a single most similar pair per family

# Output files ---------------------------------------------------------------------------------
open VOUT, ">$outfile/$category\_Vertebrate_bootstraps_$out_suffix\.txt" or die $!;
print VOUT "Ref count	Ref percent	$category count	$category percent\n";

open FOUT, ">$outfile/$category\_Fish_bootstraps_$out_suffix\.txt" or die $!;
print FOUT "Ref count	Ref percent	$category count	$category percent\n";

open KFOUT, ">$outfile/$category\_Killifish_bootstraps_$out_suffix\.txt" or die $!;
print KFOUT "Ref count	Ref percent	$category count	$category percent\n";

open WGD2ROUT, ">$outfile/$category\_2RWGD_bootstraps_$out_suffix\.txt" or die $!;
print WGD2ROUT "Ref count	Ref percent	$category count	$category percent\n";

open WGD3ROUT, ">$outfile/$category\_3RWGD_bootstraps_$out_suffix\.txt" or die $!;
print WGD3ROUT "Ref count	Ref percent	$category count	$category percent\n";

open ALLOUT, ">$outfile/$category\_All_bootstraps_$out_suffix\.txt" or die $!;
print ALLOUT "node	ref.count	ref.p	$category.count	$category.p\n";

# Open the paralog file
my @file = <PARA>;
shift @file;

# All nodes for printing -----------------------------------------------------------------
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
'Aplocheiloidei_pachypanchax',
'Nothobranchiidae_callopanchax',
'Nothobranchiidae_australe',
'Nothobranchius',
'Nothobranchius_furzeri');

# Vertebrate nodes
my %Vertebrate_nodes = ('Opisthokonta' => '',
'Bilateria' => '',
'Chordata' => '',
'Vertebrata' => '',
'Gnathostomata' => '');

# Fish nodes
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

# Killifish nodes
my %Killifish_nodes = ('Cyprinodontiformes' => '',
'Aplocheiloidei' => '',
'Aplocheiloidei_pachypanchax' => '',
'Nothobranchiidae_callopanchax' => '',
'Nothobranchiidae_australe' => '',
'Nothobranchius' => '',
'Nothobranchius_furzeri' => '');

# 10000 bootstraps
for (1..10000){

my %AgeDiffCounts;
my %AgeCounts;

my %AllCounts;
my %DiffCounts;
my $Count2R; my $Count2RDiff;
my $Count3R; my $Count3RDiff;
my $sum2r = 0; my $sum3r = 0;

# Generate an index with 50% of the data
my $fiftyP = int((scalar @file)/2);
#print "$fiftyP\n";	
# Randomize the file lines with Shuffle and pick the top 50% for each iteration
my @shuffled = shuffle(@file);		
my @picks = @shuffled[ 0 .. $fiftyP - 1 ]; 

# Foreach 50% ranndom picks, count the paralogs in each categories
foreach (@picks){ 
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
		# All count
		$AllCounts{$line[2]} += 1;
	
		# differential count
		if ($line[18] eq $category){ # Column 19 (0-based index) is category
			$DiffCounts{$line[2]} += 1;
		}	

		# 2R count
		if ($line[19] eq '2R'){ # Column 20 is 2R-WGD

			$Count2R += 1;
		
			if ($line[18] eq $category){
				$Count2RDiff += 1;
			}
		}

		# 3R count
		if ($line[20] eq '3R'){ # Column 21 is 3R-WGD
		
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
}
# SUMs
my $allsum = 0;
foreach (keys %AllCounts){$allsum = $allsum + $AllCounts{$_};}
# SUMs differential
my $diffsum = 0;
foreach (keys %DiffCounts){$diffsum = $diffsum + $DiffCounts{$_};}

# Print the counts
foreach (@all_nodes){
	
	print ALLOUT "$_\t";
	if (exists $AllCounts{$_}){
		print ALLOUT "$AllCounts{$_}\t";
		print ALLOUT eval (($AllCounts{$_}*100)/$allsum), "\t";
	}
	else {print ALLOUT "0\t0\t";}

	if (exists $DiffCounts{$_}){
		print ALLOUT "$DiffCounts{$_}\t";
		print ALLOUT eval (($DiffCounts{$_}*100)/$diffsum), "\n";
	}
	else {print ALLOUT "0\t0\n";}
		
}

# Print vertebrate bootstrap results
print VOUT $AgeCounts{'Vertebrates'}."\t";
print VOUT eval (($AgeCounts{'Vertebrates'}*100)/$allsum), "\t";
	
print VOUT $AgeDiffCounts{'Vertebrates'}."\t";
print VOUT eval (($AgeDiffCounts{'Vertebrates'}*100)/$diffsum), "\n";

# Print fish bootstrap results
print FOUT $AgeCounts{'Fish'}."\t";
print FOUT eval (($AgeCounts{'Fish'}*100)/$allsum), "\t";
	
print FOUT $AgeDiffCounts{'Fish'}."\t";
print FOUT eval (($AgeDiffCounts{'Fish'}*100)/$diffsum), "\n";

# Print killi bootstrap results
print KFOUT $AgeCounts{'Killifish'}."\t";
print KFOUT eval (($AgeCounts{'Killifish'}*100)/$allsum), "\t";
	
print KFOUT $AgeDiffCounts{'Killifish'}."\t";
print KFOUT eval (($AgeDiffCounts{'Killifish'}*100)/$diffsum), "\n";

# Print 2R bootstrap results
print WGD2ROUT "$Count2R\t";
print WGD2ROUT eval (($Count2R*100)/$allsum), "\t";
print WGD2ROUT "$Count2RDiff\t";
print WGD2ROUT eval (($Count2RDiff*100)/$diffsum), "\n";

# Print 3R bootstrap results
print WGD3ROUT "$Count3R\t";
print WGD3ROUT eval (($Count3R*100)/$allsum), "\t";
print WGD3ROUT "$Count3RDiff\t";
print WGD3ROUT eval (($Count3RDiff*100)/$diffsum), "\n";

}
print `date`;
