# Get the overlap between duplicateion file and Adam's Ka, Ks file
# There are warnings for missing Ka Ks values which is fine
use strict;
use warnings;

open OUT, ">nfur_paralogs_KaKs_31spp.txt" or die $!;

# All duplicates. File that had Ka, Ks etc 
open DUP, '../../6_paralog_SIdNdS/All_paraogs_with_SIdNdS_20200819.txt' or die $!;
my @dup = <DUP>;
my $head1 = shift @dup;
my @head1 = split "\t", $head1;
map {$_=~s/\n//g} @head1;
@head1 = @head1[6..11];
#print scalar @dup, "\n";

my %KaKs;

foreach (@dup){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;

	$KaKs{$line[1]}{$line[0]} = [@line[6..11]];
	$KaKs{$line[0]}{$line[1]} = [@line[6..11]];
	
}
#print scalar keys %KaKs, "\n";

open PARA, 'nothobranchius_furzeri_paralogs_with_counts.txt' or die $!;
my @para = <PARA>;
my $head2 = shift @para;
my @head2 = split "\t", $head2;
map {$_=~s/\n//g} @head2;
#print scalar @para;

print OUT join ("\t", @head2), "\t";
print OUT "Identitiy	Similarity	Gaps	omega	dN	dS\n";

foreach (@para){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if (exists $KaKs{$line[0]}{$line[1]}) {
	
		print OUT join ("\t", @line), "\t"; 
		print OUT join ("\t", @{$KaKs{$line[0]}{$line[1]}}), "\n";
	}
	else { # Ideally nothing should be here
		print "$line[0]\t$line[1]\n"
	}
}
print `date`;

