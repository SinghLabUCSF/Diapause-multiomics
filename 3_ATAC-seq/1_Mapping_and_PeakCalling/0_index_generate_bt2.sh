#!/bin/bash -l

# PI sunet Id or prject id to charge
#SBATCH --account=abrunet1

# Set job time to 1 day.
#SBATCH --time=1-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="drer_bt_ATAC"

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=1

# One CPU/core per task
#SBATCH --cpus-per-task=4

# 4GB of RAM
#SBATCH --mem=100G

# Who to send mail to.
#SBATCH --mail-user=param@stanford.edu

# What type of mail to send: BEGIN,END,FAIL,NONE
#SBATCH --mail-type=FAIL

# Print the stdout and stderr (they are in the same .out file) in this directory
#SBATCH --output=./log/slurm-%j.out

module load bowtie2/2.4.1
cd /labs/abrunet1/Param/ATACseq

# nfur NCBI
#bowtie2-build /labs/abrunet1/Param/RNASeq_BowTie2/genome_seq/NCBI_fasta_with_modified_header_concatenated.fa genome_index/nfur

# Aaul Dario
#bowtie2-build /labs/abrunet1/Param/RNASeq_STAR/genome_seq/aaul/AAU.v1.2.fa genome_index/aaul

# Alim
#bowtie2-build genome_seq/52670_ref_Austrofundulus_limnaeus-1.0_chrUn_renamed.fa genome_index/alim/alim

# Medaka release 100 
#bowtie2-build /labs/abrunet1/Param/RNASeq_STAR/genome_seq/olati/Oryzias_latipes.ASM223467v1.dna_sm.toplevel.fa  genome_index/olat/olat

# Zebrafish release 100
bowtie2-build /labs/abrunet1/Param/RNASeq_STAR/genome_seq/drer/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa  genome_index/drer/drer


