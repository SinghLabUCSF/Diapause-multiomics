# To sum the concentrations for each class of lipids 
# There are some warnings at line 68 that can be ignored
use strict;
#use warnings;
use List::MoreUtils 'pairwise';

# Open out file
open OUT, ">output/embryo_lipids_classes.txt" or die $!;

# Read lipid file
open LIPID, "output/embryo_lipids_normalized_Phospholipids.txt" or die $!;
my @lip = <LIPID>;
my $head = shift @lip;
my @head = split "\t", $head;
print OUT join ("\t", 'Class', @head[5..36]), "\n";

my %LipidsClasses;

# For each lipid 
foreach (@lip){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	my $class = $line[42]; # fatty acid class
	$line[41] =~s/\(|\)|\s+//g; # replace brackets from fatty acid column
	my @fa = split '/', $line[41];
	print "$line[41]\t$class\t@fa\t";
	
	# All class specific lipids
	push @{$LipidsClasses{'LIPID_CLASS'}{$class}}, \@line;

	# This part will sort the lipids into mufa, pufa or sfa for each FA and have all the intensities in values that I will sum later
	foreach (@fa){
			
		# SFA, MUFA, PUFA
		if ($_=~/.+\:0$/){
			print "sfa ";
			push @{$LipidsClasses{$class}{'SFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'SFA'}}, \@line;
		}
		if ($_=~/.+\:1$/){
			print "mufa ";
			push @{$LipidsClasses{$class}{'MUFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'MUFA'}}, \@line;
		}
		if ($_=~/.+\:([2-9])/){
			print "pufa_$1 ";
			push @{$LipidsClasses{$class}{"PUFA_$1"}}, \@line;
			push @{$LipidsClasses{$class}{"PUFA"}}, \@line;
			push @{$LipidsClasses{'ALL'}{'PUFA'}}, \@line;
		}
		
		# Ether or other linked lipids
		if ($_=~/e/){
			print "EL ";
			push @{$LipidsClasses{$class}{'EtherLinked'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'EtherLinked'}}, \@line;
		}
		if ($_=~/t/){
			print "T ";
			push @{$LipidsClasses{$class}{'TLinked'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'TLinked'}}, \@line;
		}
		if ($_=~/d/){
			print "D ";
			push @{$LipidsClasses{$class}{'DLinked'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'DLinked'}}, \@line;
		}
		
		# Compute chain length
		$_=~/(\d+)\:.+/;
		my $chain = $1;
		if ($chain <= 5){
			#print "$line[41]\t$class\t@fa\t$chain\n ";
			print "ShortChain ";
			push @{$LipidsClasses{$class}{'ShortChainFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'ShortChainFA'}}, \@line;
		}
		if ($chain >= 6 && $chain <= 12){
			print "MediumChain ";
			#print "$line[41]\t$class\t@fa\t$chain\n";
			push @{$LipidsClasses{$class}{'MediumChainFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'MediumChainFA'}}, \@line;
		}		
		if ($chain >= 13 && $chain <= 21){
			#print "$line[41]\t$class\t@fa\t$chain\n";
			print "LongChain ";
			push @{$LipidsClasses{$class}{'LongChainFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'LongChainFA'}}, \@line;
		}		
		if ($chain > 21){
			print "VeryLongChain ";
			#print "$line[41]\t$class\t@fa\t$chain\n";
			push @{$LipidsClasses{$class}{'VeryLongChainFA'}}, \@line;
			push @{$LipidsClasses{'ALL'}{'VeryLongChainFA'}}, \@line;
		}	
	}
	print "\n";
}

foreach my $class (keys %LipidsClasses){ # foreach class
	
	foreach my $fa (sort keys %{$LipidsClasses{$class}}){ # foreach lipid
		print OUT "$class\_$fa\t";
		#print scalar @{$LipidsClasses{$class}{$fa}}, "\n";
		my @lipid_lines = @{$LipidsClasses{$class}{$fa}};
		
		# Sum all lipids in this class
		my @sum; # The global sum
		foreach (@lipid_lines){
			
			my @line = @{$_}; # Only taking the samples
			my @samples = @line[5..36];
			#print "@samples\n";
			@sum = pairwise { $a + $b } @sum, @samples; # Add to the global sum. *** There are warnings here because some of the columns are text, which can be ignored
		}
		print OUT join ("\t", @sum), "\n";
	}
}

print `date`;

