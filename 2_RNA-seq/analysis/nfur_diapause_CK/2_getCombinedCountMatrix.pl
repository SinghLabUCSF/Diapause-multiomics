# Make a Exon TPM supermatrix for nfur diapause samples. Just bring counts for library and combine them into one
# This can be done with Subread etc. too, but I need to play with gene length etc. so I am dong it myself
use strict;
use warnings;

## -------------------------------------------------------------------------------------------------------------------
## For nfur_diapause_CK with NCBI genome 
my $tpmfile = '../../feature_counts';
my $outct = '1_Raw_counts/Counts_nfur_diapause_CK.csv';
my $outtpm = '1_Raw_counts/TPM_nfur_diapause_CK.csv';
my $outlen = '1_Raw_counts/GeneLengths_nfur_diapause_CK.csv';
## -------------------------------------------------------------------------------------------------------------------
## For nfur_diapause_CK with DRV genome. For more details see the GTF file source in the directory where DRV genomes are placed
# with UTR gtf file includes UTRs, so the exon TPM will exclude UTRs. This is different than NCBI GTF
#my $tpmfile = '../tpm/nfur2.0_DRV/nfur_diapause_CK_withUTRs';
#my $outfile = 'nfur2.0_DRV_diapause_CK/TPM_nfur_diapause_CK_UniquelyMapped_withUTRs.csv';
# onlyEXON GTF file, I renamed UTRs as exons, so Exon TPM will include UTRs also. Thi is similar to NCBI.
#my $tpmfile = '../tpm/nfur2.0_DRV/nfur_diapause_CK_onlyExons';
#my $outfile = 'nfur2.0_DRV_diapause_CK/TPM_nfur_diapause_CK_UniquelyMapped_onlyExons.csv';

# experimental design and libraries
my %nfurSample = ("Lib1" => "PreD1",
"Lib2" => "PreD2",
"Lib3" => "PreD3",
"Lib11" => "D3d1",
"Lib12" => "D3d2",
"Lib24" => "D3d3",
"Lib14" => "D6d1",
"Lib15" => "D6d2",
"Lib16" => "D6d3",
"Lib18" => "D1m1",
"Lib19" => "D1m2",
"Lib26" => "D1m3",
"Lib21" => "NonD1",
"Lib22" => "NonD2",
"Lib23" => "NonD3");

# Order of libraries to print
my @libOrder = ("Lib1", "Lib2", "Lib3", "Lib11", "Lib12", "Lib24", "Lib14", "Lib15", "Lib16", "Lib18", "Lib19", "Lib26", "Lib21", "Lib22", "Lib23");

my %SuperMatrix;
my %Idlist;
foreach my $lib (@libOrder){
	
	print "$tpmfile/$lib\_counts.csv\t";
	open FH, "$tpmfile/$lib\_counts.csv" or die $!;
	my @file = <FH>;
	close (FH);	
	shift @file;
	print scalar @file, "\n";
	
	foreach (@file){

		my @line = split ",", $_;
		map {$_=~s/\n//g} @line;	
		die "Lengths differ at $line[0] $lib" if ((exists $Idlist{$line[0]}) && $Idlist{$line[0]} != $line[1]);
		$Idlist{$line[0]} = $line[1]; # global id list common to all samples and their lengths
		$SuperMatrix{$line[0]}{$lib} = [$line[2], $line[3]]; # Counts and TPM
	} 
}

print "Total gene: ", scalar keys %Idlist, "\n";

# print to a file
open CT, ">$outct" or die $!;
open TPM, ">$outtpm" or die $!;
open LEN, ">$outlen" or die $!;

# Print header
foreach (@libOrder){
	if ($_ eq "Lib23"){ # If last library print without comma
		print CT "\"$nfurSample{$_}\""; print TPM "\"$nfurSample{$_}\""; print LEN "\"$nfurSample{$_}\"";
	}
	else {
		print CT "\"$nfurSample{$_}\","; print TPM "\"$nfurSample{$_}\","; print LEN "\"$nfurSample{$_}\",";
	}
}
print CT "\n"; print TPM "\n"; print LEN "\n";

foreach my $gene (keys %Idlist){
	
	print CT "\"$gene\","; print LEN "\"$gene\",";	print TPM "\"$gene\",";
	
	foreach my $lib (@libOrder){
		
		if ($lib eq "Lib23"){ # If last library print without comma
			if (exists $SuperMatrix{$gene}{$lib}){
				print LEN "$Idlist{$gene}";
				print CT "$SuperMatrix{$gene}{$lib}[0]";
				print TPM "$SuperMatrix{$gene}{$lib}[1]";

			}
			else {
				print LEN "$Idlist{$gene},";
				print CT "0,";
				print TPM "0,";			
			}
		}
		else {
			if (exists $SuperMatrix{$gene}{$lib}){
				print LEN "$Idlist{$gene},";
				print CT "$SuperMatrix{$gene}{$lib}[0],";
				print TPM "$SuperMatrix{$gene}{$lib}[1],";

			}
			else {
				print LEN "$Idlist{$gene},";
				print CT "0,";
				print TPM "0,";
			}
		}
	}
print CT "\n"; print TPM "\n"; print LEN "\n";
}
close (CT); close (TPM); close (LEN);
