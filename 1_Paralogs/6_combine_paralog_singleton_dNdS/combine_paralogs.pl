use strict;
use warnings;

my %AllParalogs;
open LARGE, '../5_orthofinder_Results_Aug04_71spp/analysis/nothobranchius_furzeri_paralogs.txt' or die $!;
my @spp71 = <LARGE>;
shift @spp71;
my %Spp71;

foreach (@spp71){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp71{$line[0]}{$line[1]} = $line[2];
	$Spp71{$line[1]}{$line[0]} = $line[2];
	
	$AllParalogs{$line[0]}{$line[1]} = '';
	$AllParalogs{$line[1]}{$line[0]} = '';	
}


open MEDIUM, '../5_orthofinder_Results_Aug11_31spp/analysis/nothobranchius_furzeri_paralogs.txt' or die $!;
my @spp31 = <MEDIUM>;
shift @spp31;
my %Spp31;

foreach (@spp31){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp31{$line[0]}{$line[1]} = $line[2];
	$Spp31{$line[1]}{$line[0]} = $line[2];
	
	$AllParalogs{$line[0]}{$line[1]} = '';
	$AllParalogs{$line[1]}{$line[0]} = '';	
}


open SMALL, '../5_orthofinder_Results_Aug12_13spp/analysis/nothobranchius_furzeri_paralogs.txt' or die $!;
my @spp13 = <SMALL>;
shift @spp13;
my %Spp13;

foreach (@spp13){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Spp13{$line[0]}{$line[1]} = $line[2];
	$Spp13{$line[1]}{$line[0]} = $line[2];
	
	$AllParalogs{$line[0]}{$line[1]} = '';
	$AllParalogs{$line[1]}{$line[0]} = '';	
}


#open ENS, '/Volumes/Mybook_3/Ohnologs/Synteny_All_2016_03_09/OhnologsForNfur/2_Paralogs/2_Nfur_paralog_annotation/8_final_nfur_duplicates/Final_nfurzeri_paralogs_new.txt' or die $!;
open ENS, '../5_Ensembl_duplicates/Final_nfurzeri_paralogs_new.txt' or die $!;
my @ens = <ENS>;
shift @ens;
my %Ens;

foreach (@ens){
	
	my @line= split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Ens{$line[0]}{$line[1]} = $line[2];
	$Ens{$line[1]}{$line[0]} = $line[2];
	
	$AllParalogs{$line[0]}{$line[1]} = '';
	$AllParalogs{$line[1]}{$line[0]} = '';
}

open OUT, ">All_paraogs_20200820.txt" or die $!;
my %checked;
print OUT "Id1	Id2	71spp	31spp	13spp	Ensembl\n";
foreach my $id1 (keys %AllParalogs){

	foreach my $id2 (keys %{$AllParalogs{$id1}}){
		
		if (not exists $checked{$id1}{$id2}){
			
			$checked{$id1}{$id2} = '';
			$checked{$id2}{$id1} = '';
			
			print OUT "$id1	$id2\t";
			if (exists $Spp71{$id1}{$id2}){
				print OUT "$Spp71{$id1}{$id2}\t";
			}
			else {print OUT "\t";}

			if (exists $Spp31{$id1}{$id2}){
				print OUT "$Spp31{$id1}{$id2}\t";
			}
			else {print OUT "\t";}

			if (exists $Spp13{$id1}{$id2}){
				print OUT "$Spp13{$id1}{$id2}\t";
			}
			else {print OUT "\t";}

			if (exists $Ens{$id1}{$id2}){
				print OUT "$Ens{$id1}{$id2}\n";
			}
			else {print OUT "\n";}		
		}
	}
}

print `date`;