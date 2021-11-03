use strict;
use warnings;

# Variables 
my $organism = 'nothobranchius_furzeri'; # The name of organism as it appears in the duplicates.tsv file
my $pattern = 'nothobranchius_furzeri_NFU\|'; # pattern for gene id as it appears in the duplicates.tsv file
my $nodes = 'nfur_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
my $outfile = 'nothobranchius_furzeri_paralogs.txt'; # Outfile with all paralogs
my $countsfile = 'nothobranchius_furzeri_paralog_counts.txt'; # Outfile with paralog counts are each nodes
my $no_paralogs = 'nothobranchius_furzeri_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

# Open outfiles
open OUT, ">$outfile" or die $!;
print OUT "gene1	gene2	duplication_node	Orthogroup	Species Tree Node	Gene Tree Node	Support	Type\n";

open OUTNO, ">$no_paralogs" or die $!;
print OUTNO "gene1	gene2	duplication_node	Orthogroup	Species Tree Node	Gene Tree Node	Support	Type	Reason for exclusion\n";

open CT, ">$countsfile" or die $!;
print  CT "duplication_node	name	genes\n";

# Read the nodes file and make a hash
open FH, $nodes or die $!;
my @nodes = <FH>;
shift @nodes;
my @node_order;
my %Nfur_nodes;

foreach (@nodes){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
		
	#print "$_";
	$Nfur_nodes{$line[0]} = $line[1];
	push @node_order, $line[0]; # order or nodes to print
}


# Read paralog file
open FH2, "../Gene_Duplication_Events/Duplications.tsv" or die $!;
#open FH2, "Duplications.test.tsv" or die $!;
my %uniqueNodes;

foreach (<FH2>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	#print "$line[0]\n";
	
	if (exists $Nfur_nodes{$line[1]}){
		
		 if (($line[5] =~/$organism/g) && ($line[6] =~/$organism/g)){
		
			#print "$line[5]\n$line[6]";
			
			my @gene1 = split ',', $line[5];
			my @gene2 = split ',', $line[6];
			map {$_=~s/\n|\s+//g} @gene1;
			map {$_=~s/\n|\s+//g} @gene2;
		
			my @notho1 = grep {$_ =~/$organism/} @gene1;
			map {$_ =~s/$pattern//g} @notho1;

			my @notho2 = grep {$_ =~/$organism/} @gene2;
			map {$_ =~s/$pattern//g} @notho2;

			#print "@notho1\t@notho2\n";
			foreach my $g1 (@notho1){
				foreach my $g2 (@notho2){
					
					print OUT "$g1	$g2	$Nfur_nodes{$line[1]}	$line[0]	$line[1]	$line[2]	$line[3]	$line[4]\n";
					$uniqueNodes{$line[1]} += 1;
				}
			}		
		}
		elsif (($line[5] =~/$organism/g) && ($line[6] !~/$organism/g)){
			print OUTNO join ("\t", @line[0..4]), "\tgene2 does not have notho gene\n";
		}
		elsif (($line[5] !~/$organism/g) && ($line[6] =~/$organism/g)){
			print OUTNO join ("\t", @line[0..4]), "\tgene1 does not have notho gene\n";
		}
		else {
			print OUTNO join ("\t", @line[0..4]), "\tboth genes do not have notho gene\n";
		}

	}
	else {
		#print "$line[1]\n";
	}
}

# Print node counts in the order
my %checked;
foreach (@node_order){
	
	print CT "$_\t$Nfur_nodes{$_}\t$uniqueNodes{$_}\n";
	if (not exists $uniqueNodes{$_}){die "$_ not found in nodes\n"}
	$checked{$_} = '';	
}

foreach (keys %uniqueNodes){
	
	if (not exists $checked{$_}){
		print CT "	$_\t$Nfur_nodes{$_}\t$uniqueNodes{$_}\n";
	}
}

print `date`;
