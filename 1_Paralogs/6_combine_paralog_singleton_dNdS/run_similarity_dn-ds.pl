use strict;
use warnings;

# Open outfile
open FINALOUT, ">All_paraogs_with_SIdNdS_20200820.txt" or die $!;

# Read final paralog file
open FH, "All_paraogs_20200820.txt" or die $!;
my @para = <FH>;
my $head = shift @para;
chomp $head;
print FINALOUT "$head\tIdentitiy	Similarity	Gaps	omega	dN	dS\n";

# Read protein file
local $/ = '>';
open PROT, "longest-Protein-seqs_Nfu-ncbi100.fasta" or die $!;
my @prot = <PROT>;
shift @prot;

my %Protein;

foreach (@prot){
	
	my @lines = split "\n", $_;
	map {$_=~s/\n//g} @lines;
	if ($lines[-1] eq '>'){pop @lines}
	#print "$lines[0]\n$lines[1]\n";
	
	my @elements = split '\|', $lines[0];
	$Protein{$elements[0]} = $lines[1];
}
print "Proteins: ", scalar keys %Protein, "\n";

# Read CDS file
open CDS, "longest-CDS-seqs_Nfu-ncbi100.fasta" or die $!;
my @cds = <CDS>;
shift @cds;

my %CDS;

foreach (@cds){
	
	my @lines = split "\n", $_;
	map {$_=~s/\n//g} @lines;
	if ($lines[-1] eq '>'){pop @lines}
	#print "$lines[0]\n$lines[1]\n";
	
	my @elements = split '\|', $lines[0];
	$CDS{$elements[0]} = $lines[1];
}
print "CDS: ", scalar keys %CDS, "\n";

# Process each pair
foreach (@para){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ((not exists $Protein{$line[0]}) || (not exists $Protein{$line[1]})){die "$line[0] or $line[1] not found in protein file\n"}
	
	print FINALOUT join ("\t", @line), "\t";
	
	# generate protein seq files
	open POUT1, ">temp/$line[0]\_prot.fasta";
	print POUT1 ">$line[0]\n$Protein{$line[0]}\n";

	open POUT2, ">temp/$line[1]\_prot.fasta";
	print POUT2 ">$line[1]\n$Protein{$line[1]}\n";
	
	# Generate nucleotide file
	open NOUT, ">temp/$line[0]\_$line[1]\_nucl.fasta";
	print NOUT ">$line[0]\n$CDS{$line[0]}\n";
	print NOUT ">$line[1]\n$CDS{$line[1]}\n";

	# Align protein and convert to nucleotide
	print `needle -asequence temp/$line[0]\_prot.fasta -bsequence temp/$line[1]\_prot.fasta -gapopen 10.0 -gapextend 0.5 -outfile temp/$line[0]\_$line[1]\_paln.fasta -aformat fasta`;
	print `needle -asequence temp/$line[0]\_prot.fasta -bsequence temp/$line[1]\_prot.fasta -gapopen 10.0 -gapextend 0.5 -outfile temp/$line[0]\_$line[1]\_paln.needle`;
	print `perl pal2nal.pl temp/$line[0]\_$line[1]\_paln.fasta temp/$line[0]\_$line[1]\_nucl.fasta -output paml -nogap -nomismatch > temp/$line[0]\_$line[1]\_caln.paml`; #  This also removes mismatch and gaps

	# Write the ctl file and run yn00
	open CTL, ">temp\/$line[0]\_$line[1]\.ctl";
	print CTL "
	seqfile = temp\/$line[0]\_$line[1]\_caln.paml * sequence data file name
    outfile = temp\/$line[0]\_$line[1]\_yn.txt * main result file
    verbose = 0 * 1: detailed output (list sequences), 0: concise output
    icode = 0 * 0:universal code; 1:mammalian mt; 2-10:see below
    weighting = 0 * weighting pathways between codons (0/1)?
    commonf3x4 = 0 * use one set of codon freqs for all pairs (0/1)? ";
	
	print `yn00 temp\/$line[0]\_$line[1]\.ctl\n`;
	
	# Parse the files
	open ALN, "temp/$line[0]\_$line[1]\_paln.needle" or die $!;

	foreach (<ALN>){
		
		if ($_=~/# Identity:\s+.+ \((.+)\%\)/g){
			print FINALOUT "$1\t";
		}
		
		if ($_=~/# Similarity:\s+.+ \((.+)\%\)/g){
			print FINALOUT "$1\t";
		}
		
		if ($_=~/# Gaps:\s+.+ \((.+)\%\)/g){
			print FINALOUT "$1\t";
		}		
	}
	
	# Parse YN output
	open YN, "temp\/$line[0]\_$line[1]\_yn.txt" or die $!;	
	foreach (<YN>){
		
		 if ($_=~/ (\d+\.?\d*) (\d+\.?\d*) \+\- \d+\.?\d*\s+(\d+\.?\d*) \+\- \d+\.?\d*\n/g){
		 	print FINALOUT "$1	$2	$3";
		 }
	}
 	print FINALOUT "\n";
				
	#exit;

}
