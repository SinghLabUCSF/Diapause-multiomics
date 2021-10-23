# This combines the information from both master files in EDGER and DESEQ2 into a single file
use strict;
use warnings;

open OUT, ">MasterPeakFile_DEseq2-EdgeR_Combined_All.txt" or die $!;

# Open master file and get all conditions for DEseq2
open MFDE, "../nfur_DEseq2/MasterPeakFile_All_DEseq2.txt" or die $!;
my @mfile = <MFDE>;
shift @mfile;

my %DEseq2;
foreach (@mfile){
	
	my @line = split "\t", $_;
	map {$_=~s/\t//g} @line;
	
	$DEseq2{$line[3]} = \@line;
}


# Open master file and get the desired conditions for EDGE-R
open MF, "../nfur_EdgeR/MasterPeakFile_All_EDGER.txt" or die $!;
my @masterfile = <MF>;
my $head = shift @masterfile;
print OUT $head;

foreach (@masterfile){
	
	my @lineOS = split "\t", $_;
	map {$_=~s/\n//g} @lineOS;
	
	my @deseq2OS = @{$DEseq2{$lineOS[3]}};
#    print "$lineOS[21]\t$lineOS[24]\t$lineOS[33]\t$lineOS[36]\t$lineOS[42]\t${$DEseq2{$lineOS[3]}}[42]\n";
#	exit;

	print OUT join ("\t", @lineOS[0..20]), "\t"; # Print first 20 columns (0 based index)
	
	Compare(21, \@lineOS, \@deseq2OS); # This is to check the column number has some value i.e. it is a DE peak or not i.e. not DE peak
	Compare(24, \@lineOS, \@deseq2OS);
	Compare(27, \@lineOS, \@deseq2OS);
	Compare(30, \@lineOS, \@deseq2OS);
	Compare(33, \@lineOS, \@deseq2OS);
	Compare(36, \@lineOS, \@deseq2OS);
	Compare(39, \@lineOS, \@deseq2OS);
	Compare(42, \@lineOS, \@deseq2OS);

	print OUT join ("\t", @lineOS[45..55]), "\n";

}
print `date`;

sub Compare {
	
	my $i = shift;

	my $ref = shift;
	my @line = @$ref;

	my $ref2 = shift;
	my @deseq2 = @$ref2;
	
	my $j = $i+2;
	
	if ($line[$i] eq '' && $deseq2[$i] eq ''){ # No DE with anyone - print EdgeR
		print OUT join ("\t", @line[$i..$j]), "\t";		
	}
	elsif ($line[$i] ne '' && $deseq2[$i] eq ''){ # DE with EdgeR; No DE with Deseq2  - print EdgeR values
		print OUT join ("\t", @line[$i..$j]), "\t";
	}
	elsif ($line[$i] eq '' && $deseq2[$i] ne ''){ 	# ************ No DE with EdgeR; DE with Deseq2  - print Deseq2 values
		print OUT join ("\t", @deseq2[$i..$j]), "\t";
	} 
	elsif ($line[$i] ne '' && $deseq2[$i] ne ''){ 	# DE with EdgeR; DE with Deseq2
		
		if ($line[$i] eq $deseq2[$i]){ # Check if the direction of change is the same - print EdgeR values
			print OUT join ("\t", @line[$i..$j]), "\t";
		}
		else { # If direction is different - die and print what's wrong.
			die "Different direction: @line[0..3]\n";
		}
	} 
	else {
		die "Not in any category: @line[0..3]\n";
	}

}
