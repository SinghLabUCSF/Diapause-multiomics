use strict;
use warnings;
use List::Util 'shuffle';

my $category = 'NeoF';

# Read the paralog file
open PARA, 'alimnaeus_DuplicatesClassification_All_le20_Reordered.txt' or die $!;

open VOUT, ">$category\_Vertebrate_bootstraps_le20.txt" or die $!;
print VOUT "Ref count	Ref percent	$category count	$category percent\n";

open FOUT, ">$category\_Fish_bootstraps_le20.txt" or die $!;
print FOUT "Ref count	Ref percent	$category count	$category percent\n";

open KFOUT, ">$category\_Killifish_bootstraps_le20.txt" or die $!;
print KFOUT "Ref count	Ref percent	$category count	$category percent\n";

=cut
open WGD2ROUT, ">$category\_2RWGD_bootstraps_le20.txt" or die $!;
print WGD2ROUT "Ref count	Ref percent	$category count	$category percent\n";

open WGD3ROUT, ">$category\_3RWGD_bootstraps_le20.txt" or die $!;
print WGD3ROUT "Ref count	Ref percent	$category count	$category percent\n";
=cut

open ALLOUT, ">$category\_All_bootstraps_le20.txt" or die $!;
print ALLOUT "node	ref.count	ref.p	$category.count	$category.p\n";

# Open the paralog file
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

# 10000 bootstraps
for (1..10000){

my %AgeDiffCounts;
my %AgeCounts;

my %AllCounts;
my %DiffCounts;
#my $Count2R; my $Count2RDiff;
#my $Count3R; my $Count3RDiff;
#my $sum2r = 0; my $sum3r = 0;

# Shuffle and take random 50%
my $fiftyP = int((scalar @file)/2);
#print "$fiftyP\n";	
my @shuffled = shuffle(@file);		
my @picks = @shuffled[ 0 .. $fiftyP - 1 ]; 

	
foreach (@picks){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
		# All count
		$AllCounts{$line[2]} += 1;
	
		# differential count
		if ($line[11] eq $category){
			$DiffCounts{$line[2]} += 1;
		}	

=cut
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
=cut		
		# Age counts
		if (exists $Vertebrate_nodes{$line[2]}){
			$AgeCounts{'Vertebrates'} += 1;
			if ($line[11] eq $category){$AgeDiffCounts{'Vertebrates'} += 1;}
		}
		if (exists $Fish_nodes{$line[2]}){
			$AgeCounts{'Fish'} += 1;
			if ($line[11] eq $category){$AgeDiffCounts{'Fish'} += 1;}
		}
		if (exists $Killifish_nodes{$line[2]}){
			$AgeCounts{'Killifish'} += 1;
			if ($line[11] eq $category){$AgeDiffCounts{'Killifish'} += 1;}
		}
}

my $allsum = 0;
foreach (keys %AllCounts){$allsum = $allsum + $AllCounts{$_};}

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

=cut
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
=cut
}


print `date`;
