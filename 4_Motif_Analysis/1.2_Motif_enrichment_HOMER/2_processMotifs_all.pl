# This is to change the motif names to gene names that make sense using clutering and motif siimlarity analysis based on TomTom
# There may be some warnings for a few directories that do not have desired files, which can be ignored.
use strict;
use warnings;

# Read the cluster Id file for motifs - from tomtom
open CL, "master_motif_similarity_tomtom/MotifClusters_filtered.txt" or die $!;
my @clust = <CL>;
shift @clust;

my %Clusters; my %GeneNames;

foreach (@clust){
	
	my @line = split "\t", $_;
	map {$_=~s/,/ /g} @line;
	#print "$line[0]\t$line[2]\t$line[4]\n";
	$Clusters{"\"$line[2]\""} = $line[0];
	$GeneNames{"\"$line[2]\""} = $line[4];
	
}

# This was run one by one for different motifs, check it for at leats one condition to see if it makes sense
# Basically I just get the motif Id from HOMER and replace with a proper curated gene name
foreach (<Motif_enrichment_results/Nfur_Master_Up/*>){
#foreach (<Motif_enrichment_results/Nfur_Master_Up/Nfur_UP_01/Either/*>){
#foreach (<Motif_enrichment_results/GenomeAlignment_Independent_Motifs/*>){
#foreach (<Motif_enrichment_results/Combined_DE_Singletons/*>){
#foreach (<Motif_enrichment_results/Positive_selection/*>){

	my $motifFile = "$_/knownResults.txt";
	my $outfile = "$_/knownResults_processed.csv";
	
	if (!-e $motifFile){		
		print "File not found: ********* Motif file: $motifFile\nOutfile: $outfile\n";
	}
	else {
		
		print "Motif file: $motifFile\nOutfile: $outfile\n\n";
		
		open OUT, ">$outfile" or die $!;

		open MOTIF, "$motifFile" or die $!;
		my @motifs = <MOTIF>;
		my $head = shift @motifs;
		my @head = split "\t", $head;
		map {$_ =~s/\n//g} @head;
		map {$_= "\"$_\""} @head;

		print OUT (join "\,", @head), ',';
		print OUT "\"FC\",\"log FC\",\"Significance(-logpvalXFC)\",\"clusterId\",\"finalGene\"\n";

		# Get the total number of motifs


		foreach (@motifs){
	
			my @line = split "\t", $_;
			#print "$line[0]\n";
			map {$_=~s/\n//g} @line;
			map {$_=~s/,/ /g} @line;
			map {$_=~s/\%//g} @line;
	
			# Compute the FC and enrichment 
			my $FC = 0;
			$FC = $line[6]/$line[8] if ($line[8] > 0);

			my $logFC = 0;
			$logFC = log10($FC) if ($FC > 0);

			my $score = 0;
			if ($line[2] == "1e-401"){$score = 401 * $logFC} # In Nfur master up the file the p-value is so small that log cannot be computed. So I am just doing it manually here.
			if ($line[2] == "1e-349"){$score = 349 * $logFC} # In Nfur master up the file the p-value is so small that log cannot be computed. So I am just doing it manually here.
			else {$score = -log10($line[2]) * $logFC};

			map {$_= "\"$_\""} @line;

			print OUT join ("\,", @line), "\,";
			print OUT "\"$FC\",\"$logFC\",\"$score\",";
			# get motif cluster
			if (exists $Clusters{$line[0]}){
				print OUT "\"$Clusters{$line[0]}\",";
			}
			else { 
				print  OUT "NA,";
			}
			
			# Print final gene name
			if (exists $GeneNames{$line[0]}){print OUT "\"$GeneNames{$line[0]}\"\n";}
			else {print "$line[0] not found in gene name\n"}
		}			
	}
#exit;
}


print `date`;


sub log10  
{ 
    my $n = shift; 
    #print "$n\n";
      
    # using pre-defined log function 
    return log($n) / log(10); 
} 

