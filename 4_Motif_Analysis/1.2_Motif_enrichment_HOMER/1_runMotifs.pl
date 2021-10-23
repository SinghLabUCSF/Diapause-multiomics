# Run motif enrichment for all lists from Adam 
use strict;
use warnings;

my %Genome = (
"Aaus" => "aaul",
"Alim" => "alim",
"Astr" => "aaul",
"Drer" => "drernew",
"Nfur" => "nfurz",
"Olat" => "olatnew");

#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_02/ConservedPeaks/REQ_1K_2O/*.bed>){  # Conserved peaks in 1 other killifish and 2 outgroups
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_02/Based_on_Alim_PCA/*.bed>){ # Motifs based on Alim PCA
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_01/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_01/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_01/Peak/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_25/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_25/Peak/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_50/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Nfur_Master_Up/Nfur_UP_50/Peak/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/lowCV_NoDE_GENES/lowCV_NoDE_GENES_25/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/lowCV_NoDE_PEAKS/lowCV_NoDE_PEAKS_01/Either/*.bed>){
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/lowCV_NoDE_PEAKS/lowCV_NoDE_PEAKS_25/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Singletons_Dia/Singletons_Dia_01/Either/*.bed>){ 
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Dev_Up_Strict/Dev_DE_01_Strict/Either/*.bed>){
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Dev_Up_Relaxed/Dev_DE_01_Relaxed/Either/*.bed>){

#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Singletons_DE_Dia/Singletons_DE_Dia_01/Either/*.bed>){
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Singletons_DE_No_Dia/Singletons_DE_No_Dia_01/*.bed>){
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Combined_DE_Singletons/Combined_Singletons_DE_01/Either/*.bed>){ # DE peaks at Singletons nfur
foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2020_03/Combined_DE_Singletons/Combined_Singletons_DE_01/Either/*.bed>){ # DE peaks at Singletons nfur

#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2021_07/Orthologs/Olat/*.bed>){ # Promoter motifs Olat
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2021_07/Orthologs/Drer/*.bed>){ # Promoter motifs Drer
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2021_07/Orthologs/Aaus/*.bed>){ # Promoter motifs Aaul
#foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2021_07/Orthologs/Alim/*.bed>){ # Promoter motifs Alim
foreach (</Volumes/Mybook_3/ATAC_Seq_Killifish/Adam_results/Bed_files/2021_07/Orthologs/Astr/*.bed>){ # Promoter motifs Astr

#	print "$_\t";

	# Get output file name
	$_ =~/\/Volumes.+\/(.+)\.bed/g;
	my $outfile = $1;
	#print "$outfile\t";

	# Get organism name for genome input
	my @org = split '_', $outfile;
	#print "$org[0]\t$org[1]\n";
	
	my $genomeName;
	if  (exists $Genome{$org[0]}){
		$genomeName = $Genome{$org[0]};
	}
	else {
		$genomeName = $Genome{$org[1]};
	}
	print "$genomeName\n";
	
#	if (($outfile =~/Neo1/) || ($outfile =~/Neo2/)){
			
		print "findMotifsGenome.pl $_ $genomeName Motif_Bed_files_Adam_results/GenomeAlignment_Independent_Motifs/$outfile -mset vertebrates\n";
		print `findMotifsGenome.pl $_ $genomeName Motif_Bed_files_Adam_results/GenomeAlignment_Independent_Motifs/$outfile -mset vertebrates`;
		#print "$outfile\n";
#	}
}

print `date`;
