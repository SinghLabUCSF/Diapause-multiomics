# To do BLAST for the input parameters provided
# 
#
#
use strict;
use warnings;

# ----------  INPUT  PARAMETERS ---------------
my $eValue = '1e-005';

if ((not defined $ARGV[0]) || (not defined $ARGV[1])){
	
	print "Usage: *.pl query_seq_path database_path\n\n";
	print "Please provide 2 paths: ARGV[0] => path to a query sequence file, and ARGV[1] => path to blast database\n\n";
	print "*** IMPORTANT: Both the inputs must be a path. If in the same directory, use full paths\n\n";
	die;
}

$ARGV[0] =~/.+\/(.+)\..+/g;
#print "$1\n";
my $from = $1;

$ARGV[1] =~/.*\/(.+)/g;
#print "$1\n";
my $to = $1;

#print "Sequence: $from\nDatabase: $to\n";

my $querySeqFile = $ARGV[0];
my $database = $ARGV[1];
my $bestHitOutfile = "BestHits_$from\-to-$to\.txt";


# Read sequence fasta file
local $/ = '>';
open FH, $querySeqFile or die $!;
my @sequenceFile = <FH>;
close (FH);
shift (@sequenceFile);

# open best hit and all hit output file
open OUT, ">$bestHitOutfile" or die $!;
#open OUT2, ">$allHitsOutfile" or die $!;

# for each sequence generate a temporary sequence file to perform BLAST and parse the best hit
my $count = 1;
foreach (@sequenceFile){
	
	my @lines = split "\n", $_;
	if ($lines[-1] eq '>'){pop @lines}
	
	# generate temporary sequence file having one sequence
	my $tempseqfile = "$from\-$to\-tempSeqFile.txt";
	my $tempout = "$from\-$to\-tempout.txt";

	open TEMP, ">$tempseqfile" or die $!;
	print TEMP ">";
	print TEMP join "\n", @lines;
	close (TEMP);
	
	# perfoem BLAST against all human proteins db (evalue 1e-007)
	`blastp -query $tempseqfile -task blastp -db $database -out $tempout -evalue $eValue -outfmt "6 qseqid sseqid evalue bitscore score length"`;
	
	# parse the outputs
	parseBestHit($tempout);
	
	print 'Processed '.$count.' of '.scalar (@sequenceFile)."\n";
	$count++;
	
}

# parse BLAST hit

sub parseBestHit {
	
	my $tmpout = shift;
	# get all the hits in an array
	local $/ = "\n";
	open FH2, "$tmpout" or die $!;
	my @out = <FH2>;
	
	# print best hit in output file if there is a best hit below given e-value
	print OUT $out[0] if (exists $out[0]);
	
	# print all hits below evalue to another file
	#print OUT2 @out;
}


