# This is to get the duplicates that are not in Ensembl but a part of the gene family clustering.
# The node names are from Ensembl and NCBI Taxonomy
# Family clustering is from Wagner et al. 2018
use strict;
use warnings;

local $| = 1;

my $ensembl_duplicates = '5_Combined_nodes/nfurzeri_consensus_dup_nodes.txt';
my $fish_families = '/Volumes/Mybook_2/Coding_genome_analysis_with_Alim/4_Protein_clusters/proteinortho_run_10-3-16_2.proteinortho';
open OUT, '>7_get_duplication_time_from_families/nfur_duplicates_from_families.txt' or die $!;
print OUT "Id1	Id2	Duplicate in other Notho	Duplicate in non-annual	Duplicate in other teleosts	Assigned node\n";

# make duplicate hash in both directions ---------------------------------------
open DUP, "$ensembl_duplicates" or die $!;
my @dupli = <DUP>;

my %Duplicates;

foreach (@dupli){
	
	my @line = split "\t", $_;
	map {$_ =~s/\s//g} @line;

	$Duplicates{$line[0]}{$line[1]} = $line[2];
	$Duplicates{$line[1]}{$line[0]} = $line[2];
}
#print scalar keys %Duplicates;

# Define hash to count duplicates in other species -----------------------------
# Other Notho fishes from transriptome 
my %Nothos = (      '5' => 'NPI',
                    '7' => 'NKO',
                    '8' => 'NKA',
                    '10' => 'NRA',
                    '12' => 'NKU');
# Non annual Ast or Alim
my %NonAnnual = (   '6' => 'AST',
					'19' => 'ALIM',
);
# Other distantly related teleost fish species                    
my %OtherTeleosts= ('3' => 'gmorhua',
                    '4' => 'drerio',
                    '9' => 'tnigroviridis',
                    '11' => 'olatipes',
                    '13' => 'gaculeatus',
                    '14' => 'trubripes',
                    '15' => 'xmaculatus',
                    '16' => 'oniloticus',
                    '17' => 'amexicanus',
                    '18' => 'pformosa',
                    '21' => 'fheteroclitus');

# Read protein ortho file and get the duplicates in Nfur ------------------------
open FAM, "$fish_families" or die $!;
my @families = <FAM>;
shift @families;

my %Marked; # To check in both direction

foreach (@families){
	
	my @line = split "\t", $_;
	map {$_ =~s/\s//g} @line;	
	
	if ($line[20] =~/\,/g){
		
		$line[20] =~s/NFU\|//g;
#		print ">$line[20]\n";
		my @dup = split ',', $line[20];
		
		for (my $i = 0; $i < $#dup; $i++){
			for (my $j = $i+1; $j <= $#dup; $j++){

#				print "\t$dup[$i]\t$dup[$j]\n";
				if ((not exists $Duplicates{$dup[$i]}{$dup[$j]}) && (not exists $Marked{$dup[$i]}{$dup[$j]})){
				
					print OUT "$dup[$i]\t$dup[$j]\t";
					$Marked{$dup[$i]}{$dup[$j]} = '';
					$Marked{$dup[$j]}{$dup[$i]} = '';
					
					# Count how many teleosts, non-annual and notho species also have duplication					
					my $teleostCount = 0; my $otherNothoCount = 0; my $nonAnnualCount = 0;
					
					# Find out how many teleosts are there
					foreach (keys %OtherTeleosts){
						if ($line[$_] =~/\,/g){$teleostCount++;}
					}
					# Find out the number of nothos except Nfur
					foreach (keys %Nothos){
						if ($line[$_] =~/\,/g){$otherNothoCount++;}
					}
					# Find out the number of non annuals
					foreach (keys %NonAnnual){
						if ($line[$_] =~/\,/g){$nonAnnualCount++;}
					}
					
					print OUT "$otherNothoCount\t$nonAnnualCount\t$teleostCount\t";
					
					# Decide node based on where the duplication happened
					# If teleosts have a duplication -> After teleosts
					# Elsif teleost don't but non-annual do -> In non annual
					# Elsif both teleost and non-annual don't but notho. do -> In notho branch
					# Else if none of them have -> Specifically in nothobranchius
					
					if ($teleostCount > 0){print OUT "Teleosts"}
					elsif ($nonAnnualCount > 0){print OUT "Aplocheiloidei"}
					elsif ($otherNothoCount > 0){print OUT "Nothobranchiadae"}
					else {print OUT "Nothobranchius furzeri"}
					print OUT "\n";
				}
			}
		}
		#if (scalar @dup == 2){
		#	if (not exists $Duplicates{$dup[0]}{$dup[1]}){print "$line[20]\n";}
		#}
	}
}

print `date`;



