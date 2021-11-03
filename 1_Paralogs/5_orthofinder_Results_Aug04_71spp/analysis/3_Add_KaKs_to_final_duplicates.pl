# Get the overlap between duplicateion file and protein evolutionary rate or Ka, Ks file.
# There are warnings for missing Ka, Ks values which is fine.
use strict;
use warnings;

open OUT, ">nfur/nfur_paralogs_KaKs_71spp.txt" or die $!;

# All duplicates combined from all analyses with Ka, Ks etc.
open DUP, '../../6_combine_paralog_singleton_dNdS/All_paraogs_with_SIdNdS_20200819.txt' or die $!;
my @dup = <DUP>;
my $head1 = shift @dup;
my @head1 = split "\t", $head1;
map {$_=~s/\n//g} @head1;
@head1 = @head1[6..11];
#print scalar @dup, "\n";

my %KaKs;

foreach (@dup){ # Read the ka, ks file and make a hash with relevant info in column 6..11
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;

	$KaKs{$line[1]}{$line[0]} = [@line[6..11]]; # Making a hash in both directions to cover the pair directinality
	$KaKs{$line[0]}{$line[1]} = [@line[6..11]];
	
}
#print scalar keys %KaKs, "\n";

# Read the paralog file in which I have to add the Ka, Ks information
open PARA, 'nfur/nothobranchius_furzeri_paralogs_with_counts.txt' or die $!;
my @para = <PARA>;
my $head2 = shift @para;
my @head2 = split "\t", $head2;
map {$_=~s/\n//g} @head2;
#print scalar @para;

print OUT join ("\t", @head2), "\t";
print OUT "Identitiy	Similarity	Gaps	omega	dN	dS\n";

foreach (@para){ # Foreach para
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if (exists $KaKs{$line[0]}{$line[1]}) { # If it exists in the Ka Ks file
	
		print OUT join ("\t", @line), "\t"; 
		print OUT join ("\t", @{$KaKs{$line[0]}{$line[1]}}), "\n"; # Print it's infomation
	}
	else { # Ideally nothing should be here because all paralogs should be in the combined file
		print "$line[0]\t$line[1]\n"
	}
}
print `date`;
