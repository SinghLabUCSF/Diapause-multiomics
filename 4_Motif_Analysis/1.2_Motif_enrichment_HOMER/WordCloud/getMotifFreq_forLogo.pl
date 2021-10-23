# Get frequency for diapause specific motifs to plot as wordart
# I will read diapause motifs, exclude development specific motifs, and count occurrances based on P-value
use strict;
use warnings;

open OUT, ">Diapause_specific_motif_freq.txt" or die $!;

# Development motifs
open DEVMOT, "Strict_All-Development-UP/knownResults_processed.csv" or die $!;
my @devmotif = <DEVMOT>;
shift @devmotif;
my %DevMotif;

foreach (@devmotif){
	
	my @line = split ',', $_;
	map {$_=~s/\"|\n//g} @line;
	
	if ($line[4] < 0.1){
		$DevMotif{$line[13]} = \@line;
	}
}
#print scalar keys %DevMotif;

# Diapause motifs
open DIAMOT, "../Motif_enrichment_results/Nfur_Master_Up/Nfur_Neo_Combine/knownResults_processed.csv" or die $!;
my @motif = <DIAMOT>;
shift @motif;
my %MotifCluster;
foreach (@motif){
	
	my @line = split ',', $_;
	map {$_=~s/\"|\n//g} @line;
	
	# Get the -log pvalue
	my $decimal = - sprintf("%.10g", $line[3]);
	
	# Print motifs based on -log P for logo if they don't exist on developmental motifs
	if (not exists $DevMotif{$line[13]}){
		
		# clean up a bit
		$line[13] =~s/\(.+\)//g;
		
		for (0..$decimal){
			
			print OUT "$line[13]\n";
#			if (exists $Family{$line[0]}){
#				print "$Family{$line[0]}\n";
#			}
#			else {print "**** Error: $line[0]\n"}
#			
			$MotifCluster{$line[12]}{$line[13]} += 1;
		}
	}
}

foreach my $clust (keys %MotifCluster){
	foreach (keys %{$MotifCluster{$clust}}){
#		print "$clust\t$_\t$MotifCluster{$clust}{$_}\n";
	}
}

print `date`