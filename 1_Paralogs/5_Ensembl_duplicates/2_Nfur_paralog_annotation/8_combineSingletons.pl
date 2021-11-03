# Combine singleton genes from all 4 fish - and check that these genes are not in the duplication file
# Output: A file with all genes that are singletons wrt at least 1 of the fish, and are not in the final duplication file

use strict;
use warnings;

# Read the final combined Paralog file
my %Paralogs;

open PARA, "8_final_nfur_duplicates\/Final_nfurzeri_paralogs_new.txt" or die $!;
my @paralogs = <PARA>;
shift @paralogs;

foreach (@paralogs){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	#print "$line[0]\t$line[1]\n";
	$Paralogs{$line[0]} = '';
	$Paralogs{$line[1]} = '';	
}
#print scalar keys %Paralogs;

# Outfile 
open OUT, '>8_final_nfur_duplicates/Final_nfurzeri_singletons.txt' or die $!; # 3R

my %Singletons;

foreach (<6_Sigletons\/nfurzeri_singletons_from_*.txt>){
	
	print "$_\n";
	$_=~/.+\/nfurzeri_singletons_from_(.+)\.txt/g;
	my $org = $1;
	#print "$org\n";
	
	open FH, "$_" or die $!;
	my @singletons = <FH>;
	shift @singletons;	
	
	foreach (@singletons){
		
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
		
	
		if ((not exists $Singletons{$line[0]}) && (not exists $Paralogs{$line[0]})){
		
			print OUT "$line[0]\n";
			$Singletons{$line[0]} = '';
		}
	}
}



