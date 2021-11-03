# Parse the OrthFInder results and generate a table of paralog and their duplication time
use strict;
use warnings;

# Run it for different organisms one by one. Check for at least one species, preferably Nfur
# Variables
my $organism = 'nothobranchius_furzeri'; # The name of organism as it appears in the duplicates.tsv file
my $pattern = 'nothobranchius_furzeri_NFU\|'; # pattern for gene id as it appears in the duplicates.tsv file
my $nodes = 'nfur/nfur_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
my $outfile = 'nfur/nothobranchius_furzeri_paralogs.txt'; # Outfile with all paralogs *** This will be used for further analysis ***
my $countsfile = 'nfur/nothobranchius_furzeri_paralog_counts.txt'; # Outfile with paralog counts are each nodes
my $no_paralogs = 'nfur/nothobranchius_furzeri_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

#my $organism = 'austrofundulus_limnaeus'; # The name of organism as it appears in the duplicates.tsv file
#my $pattern = 'austrofundulus_limnaeus_alimnaeus\|'; # pattern for gene id as it appears in the duplicates.tsv file
#my $nodes = 'alim_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
#my $outfile = 'austrofundulus_limnaeus_paralogs.txt'; # Outfile with all paralogs
#my $countsfile = 'austrofundulus_limnaeus_paralog_counts.txt'; # Outfile with paralog counts are each nodes
#my $no_paralogs = 'austrofundulus_limnaeus_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

#my $organism = 'danio_rerio'; # The name of organism as it appears in the duplicates.tsv file
#my $pattern = 'danio_rerio_'; # pattern for gene id as it appears in the duplicates.tsv file
#my $nodes = 'zebra_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
#my $outfile = 'danio_rerio_paralogs.txt'; # Outfile with all paralogs
#my $countsfile = 'danio_rerio_paralog_counts.txt'; # Outfile with paralog counts are each nodes
#my $no_paralogs = 'danio_rerio_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

#my $organism = 'oryzias_latipes'; # The name of organism as it appears in the duplicates.tsv file
#my $pattern = 'oryzias_latipes_'; # pattern for gene id as it appears in the duplicates.tsv file
#my $nodes = 'medaka_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
#my $outfile = 'oryzias_latipes_paralogs.txt'; # Outfile with all paralogs
#my $countsfile = 'oryzias_latipes_paralog_counts.txt'; # Outfile with paralog counts are each nodes
#my $no_paralogs = 'oryzias_latipes_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

#my $organism = 'aphyosemion_australe'; # The name of organism as it appears in the duplicates.tsv file
#my $pattern = 'aphyosemion_australe_'; # pattern for gene id as it appears in the duplicates.tsv file
#my $nodes = 'aaul_nodes.txt'; # Nodes leading to organism under investigation, manually extracted from species tree
#my $outfile = 'aphyosemion_australe_paralogs.txt'; # Outfile with all paralogs
#my $countsfile = 'aphyosemion_australe_paralog_counts.txt'; # Outfile with paralog counts are each nodes
#my $no_paralogs = 'aphyosemion_australe_missed.txt'; # Outfile with missed paralogs because both the genes don't match the pattern

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


# Read paralog file from orthoFInder results.
# Duplications.tsv has all the duplicates for all species in a single file in column 5 and 6 (Perl column index starts from 0)
open FH2, "../Gene_Duplication_Events/Duplications.tsv" or die $!;
#open FH2, "Duplications.test.tsv" or die $!;
my %uniqueNodes;

foreach (<FH2>){ # Foreach line
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	#print "$line[0]\n";
	
	# Read only the nodes for Nfur, to ignore other species
	if (exists $Nfur_nodes{$line[1]}){
		
		 if (($line[5] =~/$organism/g) && ($line[6] =~/$organism/g)){ # If the paralog line includes the current organism, then only process, else it is for other species
		
			#print "$line[5]\n$line[6]";
			
			my @gene1 = split ',', $line[5]; # get individual gene ids in an array
			my @gene2 = split ',', $line[6];
			map {$_=~s/\n|\s+//g} @gene1; # Remove \n etc.
			map {$_=~s/\n|\s+//g} @gene2;
			
			# Grep all the ids for the current organism into arrays
			my @notho1 = grep {$_ =~/$organism/} @gene1;
			map {$_ =~s/$pattern//g} @notho1; # Remove unwanted things

			my @notho2 = grep {$_ =~/$organism/} @gene2; 
			map {$_ =~s/$pattern//g} @notho2;

			#print "@notho1\t@notho2\n";
			# Now make paralog pairs and print with relevant info
			foreach my $g1 (@notho1){
				foreach my $g2 (@notho2){
					
					print OUT "$g1	$g2	$Nfur_nodes{$line[1]}	$line[0]	$line[1]	$line[2]	$line[3]	$line[4]\n";
					$uniqueNodes{$line[1]} += 1; # This is to get a list of unique node names
				}
			}		
		} # The things below are saniuty check for other patterns
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

# Print counts for each node in the order of their age
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
