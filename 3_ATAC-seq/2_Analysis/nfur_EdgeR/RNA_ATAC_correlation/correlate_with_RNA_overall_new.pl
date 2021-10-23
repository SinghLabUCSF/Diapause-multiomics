# Get median RNA and ATAcseq values
use strict;
use warnings;
use List::Util qw(min max);
use Statistics::Basic qw(:all);

# RNA-seq
open RNA, "Counts_DE_Combined_Final.csv" or die $!; # RNA-seq counts
my @rna = <RNA>;
shift @rna;

my %rnaDIA; my %rnaDEV;

foreach (@rna){
	
	my @line = split ',', $_;
	map {$_ =~s/\n|\"//g} @line;
	#print "$line[0]\n";
	
	$rnaDIA{$line[0]} = $line[32]; # Diapause median RNA level from all replicates and conditions
	$rnaDEV{$line[0]} = $line[31]; # Development median RNA level from all replicate and conditions
	
}


# Annotation and ATAC
open ATAC, "../consensusPeaks_annotations.txt" or die $!; # ATAC-seq values
my @atac = <ATAC>;
shift @atac;

# Get diapause values
my %atacDIA;

foreach (@atac){
	
	my @line = split "\t", $_;
	map {$_ =~s/\n//g} @line;
	#print "$line[25]\n";
	
	if ($line[19] =~ /Promoter/g){
		# Column 25 (starting from 0) is gene name
		# Since there can be multiple peaks for each gene and replicate across multiple lines. This will put them all in an array. Then I will take median or maximum.
		push @{$atacDIA{$line[25]}}, $line[8]; # dia 6 day rep 1
		push @{$atacDIA{$line[25]}}, $line[9]; # dia 6d rep 2
		push @{$atacDIA{$line[25]}}, $line[15]; # dia 1m rep 1
		push @{$atacDIA{$line[25]}}, $line[16]; # dia 1m rep 2
		push @{$atacDIA{$line[25]}}, $line[17]; # dia 1m rep 4
	}	
}
#print join "\n", @{$atacDIA{"fmr1"}},"\n";

# Get development values
my %atacDEV;

foreach (@atac){
	
	my @line = split "\t", $_;
	map {$_ =~s/\n//g} @line;
	#print "$line[24]\n";
	
	# Column 25 is gene name
	# Since there can be multiple peaks for each gene and replicate across multiple lines. This will put them all in an array. Then I will take median or maximum.
	push @{$atacDEV{$line[25]}}, $line[10]; # dev rep 1
	push @{$atacDEV{$line[25]}}, $line[11]; # dev rep 2
	push @{$atacDEV{$line[25]}}, $line[13]; # dev rep 1
	push @{$atacDEV{$line[25]}}, $line[14]; # dev rep 2
	
}
#print join "\n", @{$atacDEV{"fmr1"}},"\n";


# Now take median and max of all peaks
# for diapsuse
my %atacDIAmean; my %atacDIAmax; 

foreach my $gene (keys %atacDIA){
	my $mean = median(@{$atacDIA{$gene}});
	$mean =~s/\,//g;
	$atacDIAmean{$gene} = $mean;
	
	my $max = max(@{$atacDIA{$gene}});
	$max =~s/\,//g;
	$atacDIAmax{$gene} = $max;
}
#print $atacDIAmean{"fmr1"},"\n";
#print print $atacDIAmax{"fmr1"};

# for developemnt
my %atacDEVmean; my %atacDEVmax;

foreach my $gene (keys %atacDEV){

	my $mean = median(@{$atacDEV{$gene}});
	$mean =~s/\,//g;
	$atacDEVmean{$gene} = $mean;
	
	my $max = max(@{$atacDEV{$gene}});
	$max =~s/\,//g;
	$atacDEVmax{$gene} = $max;
}


# Now combine ATAC and RNA. It is expected that RNA will have more  genes so I start with RNA and then add ATAC data
my %RNA_ATAC_combined;

# For diapause
open OUTDIA, ">RNA-ATAC_correlation_diapause_all.txt" or die $!;
print OUTDIA "gene	rna	atac_median	atac_max\n";
foreach (keys %rnaDIA){
	
	if ((exists $atacDIAmean{$_}) && (exists $atacDIAmax{$_})){
#		if ($rnaDIA{$_}	> 10 && $atacDIAmean{$_} > 5){
			print OUTDIA "$_	$rnaDIA{$_}	$atacDIAmean{$_}	$atacDIAmax{$_}\n";
#		}
	}	
}

# For development
open OUTDEV, ">RNA-ATAC_correlation_development_all.txt" or die $!;
print OUTDEV "gene	rna	atac_median	atac_max\n";
foreach (keys %rnaDEV){
	
	if ((exists $atacDEVmean{$_}) && (exists $atacDEVmax{$_})){
		
#		if ($rnaDEV{$_}	> 10 && $atacDEVmean{$_} > 5){
			print OUTDEV "$_	$rnaDEV{$_}	$atacDEVmean{$_}	$atacDEVmax{$_}\n";
		}
#	}
}

print `date`;



