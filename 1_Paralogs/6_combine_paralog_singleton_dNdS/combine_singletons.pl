use strict;
use warnings;

my %AllSingletons;
open LARGE, '../5_orthofinder_Results_Aug04_71spp/analysis/nfur/nothobranchius_furzeri_singletons.txt' or die $!;
my @spp71 = <LARGE>;
shift @spp71;
my %Spp71;

foreach (@spp71){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp71{$line[1]} = $line[0];
	
	$AllSingletons{$line[1]} = '';
}


open MEDIUM, '../5_orthofinder_Results_Aug11_31spp/analysis/nothobranchius_furzeri_singletons.txt' or die $!;
my @spp31 = <MEDIUM>;
shift @spp31;
my %Spp31;

foreach (@spp31){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp31{$line[1]} = $line[0];
	
	$AllSingletons{$line[1]} = '';
}


open SMALL, '../5_orthofinder_Results_Aug12_13spp/analysis/nothobranchius_furzeri_singletons.txt' or die $!;
my @spp13 = <SMALL>;
shift @spp13;
my %Spp13;

foreach (@spp13){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp13{$line[1]} = $line[0];
	
	$AllSingletons{$line[1]} = '';
}


#open ENS, '/Volumes/Mybook_3/Ohnologs/Synteny_All_2016_03_09/OhnologsForNfur/2_Paralogs/2_Nfur_paralog_annotation/8_final_nfur_duplicates/Final_nfurzeri_paralogs_new.txt' or die $!;
open ENS, '../5_Ensembl_duplicates/Final_nfurzeri_singletons.txt' or die $!;
my @ens = <ENS>;
my %Ens;

foreach (@ens){
	
	$_=~s/\n//g;
	
	$Ens{$_} = '';
	
	$AllSingletons{$_} = '';
}

open OUT, ">All_singletons_20210111.txt" or die $!;
my %checked;
print OUT "Id	71spp	31spp	13spp	Ensembl\n";
foreach my $id (keys %AllSingletons){
		
	print OUT "$id\t";
	if (exists $Spp71{$id}){
		print OUT "Y\t";
	}
	else {print OUT "\t";}

	if (exists $Spp31{$id}){
		print OUT "Y\t";
	}
	else {print OUT "\t";}

	if (exists $Spp13{$id}){
		print OUT "Y\t";
	}
	else {print OUT "\t";}

	if (exists $Ens{$id}){
		print OUT "Y\n";
	}
	else {print OUT "\n";}
}

print `date`;