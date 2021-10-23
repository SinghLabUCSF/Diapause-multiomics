#!/bin/bash -l

# PI sunet Id or prject id to charge
#SBATCH --account=abrunet1

# Set job time to 1 day. One day:--time=1-00:00:00
#SBATCH --time=4-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="star"

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=1

# One CPU/core per task
#SBATCH --cpus-per-task=2

# 100GB of RAM. 100G typically
#SBATCH --mem=100G

# Who to send mail to.
#SBATCH --mail-user=param@stanford.edu

# What type of mail to send: BEGIN,END,FAIL,NONE
#SBATCH --mail-type=NONE

# Print the stdout and stderr (they are in the same .out file) in this directory
#SBATCH --output=./log/slurm-%j.out

cd /labs/abrunet1/Param/ATACseq
module load trim_galore/0.4.5
module load bowtie2/2.4.1
module load samtools/1.5
module load picard/2.22.8
module load deeptools/3.3.0
module load macs2/2.2.71
module load bedtools/2.27.1


# TRIM-GALORE: adoptor trimming =======================================================================================================================
##echo "trim_galore --fastqc --fastqc_args \"--outdir fastqc/$1\" --gzip --output_dir $6 --paired $6/$2 $6/$3"
#trim_galore --fastqc --fastqc_args "--outdir fastqc/$1" --gzip --output_dir $6 --paired $6/$2 $6/$3

# BowTie2 read alignment. LOG is wriien in QC DIRECTORY.================================================================================================
##echo "bowtie2 --very-sensitive -x $4 -1 $6/$5\_R1_001_val_1.fq.gz -2 $6/$5\_R2_001_val_2.fq.gz -S alignment_SAM/$1/$5\.sam 2>&1 | tee QC/$1/$5\_alignment.log"
# FOR: NFUR, AAUL, ALIM: bowtie2 --very-sensitive -x $4 -1 $6/$5\_R1_001_val_1.fq.gz -2 $6/$5\_R2_001_val_2.fq.gz -S alignment_SAM/$1/$5\.sam 2>&1 | tee QC/$1/$5\_alignment.log
# FOR: All other genomes bowtie2 --very-sensitive -x $4 -1 $6/$5\_1_val_1.fq.gz -2 $6/$5\_2_val_2.fq.gz -S alignment_SAM/$1/$5\.sam 2>&1 | tee QC/$1/$5\_alignment.log

#-> starting here I am commenting to run only bw files using #->
#->if [[ ("$5" =~ "Olat") || ("$5" =~ "Drer")]]; then
#->    echo "Drer or Olati"
#->    bowtie2 --very-sensitive -x $4 -1 $6/$5\_1_val_1.fq.gz -2 $6/$5\_2_val_2.fq.gz -S alignment_SAM/$1/$5\.sam 2>&1 | tee QC/$1/$5\_alignment.log
#->else
#->    echo "killifish"
#->    bowtie2 --very-sensitive -x $4 -1 $6/$5\_R1_001_val_1.fq.gz -2 $6/$5\_R2_001_val_2.fq.gz -S alignment_SAM/$1/$5\.sam 2>&1 | tee QC/$1/$5\_alignment.log
#->fi

# convert SAM file to BAM file
##echo "samtools view -S -b alignment_SAM/$1/$5\.sam > alignment_BAM/$1/$5\.bam"
#->samtools view -S -b alignment_SAM/$1/$5\.sam > alignment_BAM/$1/$5\.bam

# Sort filtered Bam based on position
# I am using input file name as temp file prefix to avoid any issues due to overlapping temp file name.
##echo "samtools sort -T alignment_BAM/$1/$5\.bam -o alignment_BAM/$1/$5\_sorted.bam alignment_BAM/$1/$5\.bam" 
#->samtools sort -T alignment_BAM/$1/$5\.bam -o alignment_BAM/$1/$5\_sorted.bam alignment_BAM/$1/$5\.bam
####echo "samtools sort -n -T alignment_BAM/$1/$5\.bam -o alignment_BAM/$1/$5\_sorted.bam alignment_BAM/$1/$5\.bam" # sort based on name for Genrich

# Index bam file
##echo "samtools index alignment_BAM/$1/$5\_sorted.bam"
#->samtools index alignment_BAM/$1/$5\_sorted.bam

# Count mito reads
#->echo "100 * $(samtools view -c alignment_BAM/$1/$5\_sorted.bam $7) / $(samtools view -c alignment_BAM/$1/$5\_sorted.bam)" | bc -l > QC/$1/$5\.mito_read_percent.txt

# Remove mito reads
##echo "samtools view -h alignment_BAM/$1/$5\_sorted.bam | grep -v $7 | samtools sort -O bam -o alignment_BAM/$1/$5\_mitoRemoved.bam"
#->samtools view -h alignment_BAM/$1/$5\_sorted.bam | grep -v $7 | samtools sort -O bam -o alignment_BAM/$1/$5\_mitoRemoved.bam

# Mark duplicates. Matrix file which has duplicate infor is put in QC folder, and BAM file in BAM folder
##echo "picard MarkDuplicates QUIET=true INPUT=alignment_BAM/$1/$5\_mitoRemoved.bam OUTPUT=alignment_BAM/$1/$5\_dupMarked.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=QC/$1/$5\_dupMarked.sorted.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$5\_tmp"
#->export _JAVA_OPTIONS='-Xms10G -Xmx15G -XX:ParallelGCThreads=1'
#->picard MarkDuplicates QUIET=true INPUT=alignment_BAM/$1/$5\_mitoRemoved.bam OUTPUT=alignment_BAM/$1/$5\_dupMarked.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=QC/$1/$5\_dupMarked.sorted.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT  TMP_DIR=$5\_tmp

# Remove duplicates
##echo "samtools view -h -b -F 1024 alignment_BAM/$1/$5\_dupMarked.bam > alignment_BAM/$1/$5\_dupRemoved.bam"
#->samtools view -h -b -F 1024 alignment_BAM/$1/$5\_dupMarked.bam > alignment_BAM/$1/$5\_dupRemoved.bam

# Remove multi-mapped reads (i.e. those with MAPQ < 20, using -q in SAMtools)
##echo "samtools view -h -q 20 alignment_BAM/$1/$5\_dupRemoved.bam > alignment_BAM/$1/$5\_mapQ20.bam"
#->samtools view -h -q 20 alignment_BAM/$1/$5\_dupRemoved.bam > alignment_BAM/$1/$5\_mapQ20.bam

# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2
##echo "samtools view -h -b -F 1804 -f 2 alignment_BAM/$1/$5\_mapQ20.bam > alignment_BAM/$1/$5\_filtered.bam"
#->samtools view -h -b -F 1804 -f 2 alignment_BAM/$1/$5\_mapQ20.bam > alignment_BAM/$1/$5\_filtered.bam

# Shift alignments using alignmentSieve in DeepTools
##echo "samtools sort -T alignment_BAM/$1/$5\_filtered.bam -o alignment_BAM/$1/$5\_filtered_sorted.bam alignment_BAM/$1/$5\_filtered.bam"
#->samtools sort -T alignment_BAM/$1/$5\_filtered.bam -o alignment_BAM/$1/$5\_filtered_sorted.bam alignment_BAM/$1/$5\_filtered.bam

##echo "samtools index alignment_BAM/$1/$5\_filtered_sorted.bam"
#->samtools index alignment_BAM/$1/$5\_filtered_sorted.bam

##echo "alignmentSieve --numberOfProcessors 4 --ATACshift --bam alignment_BAM/$1/$5\_filtered_sorted.bam -o alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam"
#->alignmentSieve --numberOfProcessors 4 --ATACshift --bam alignment_BAM/$1/$5\_filtered_sorted.bam -o alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam

# the bam file needs to be sorted and indexed again
##echo "samtools sort -O bam -o alignment_BAM/$1/$5\_shifted_FINAL.bam alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam"
#->samtools sort -O bam -o alignment_BAM/$1/$5\_shifted_FINAL.bam alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam

##echo "samtools index alignment_BAM/$1/$5\_shifted_FINAL.bam"
#->samtools index alignment_BAM/$1/$5\_shifted_FINAL.bam
##echo "rm alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam"
#->rm alignment_BAM/$1/$5\_filtered_sorted.bam\.tmp.bam

# Convet to bigwig - I tried different things and saw that bin 1 and rpkm normalization produces best looking peaks
#echo "bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL.bw --extendReads"
#bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL.bw --extendReads
#bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL_rpkm_bin10.bw --extendReads 100 --normalizeUsing RPKM --binSize 10
bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL_rpkm_bin1.bw --extendReads 100 --normalizeUsing RPKM --binSize 1
#bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL_rpgc_bin10.bw --extendReads 100 --normalizeUsing RPGC --effectiveGenomeSize $9 --binSize 10
#bamCoverage -b alignment_BAM/$1/$5\_shifted_FINAL.bam -o alignment_BAM/$1/$5\_shifted_FINAL_rpgc_bin1.bw --extendReads 100 --normalizeUsing RPGC --effectiveGenomeSize $9 --binSize 1

# Call peaks
##echo "MACS2 callpeak -f BAMPE -g $9 --keep-dup all --cutoff-analysis -n $5 -B -t alignment_BAM/$1/$5\_shifted_FINAL.bam --outdir peaks/$1  2> peaks/$1/$5\_macs2.log"
#->macs2 callpeak -f BAMPE -g $9 --keep-dup all --cutoff-analysis -n $5 -B -t alignment_BAM/$1/$5\_shifted_FINAL.bam --outdir peaks/$1  2> peaks/$1/$5\_macs2.log


## QC and heatmaps etc.
# The part below is to get the TSS coverage heatmaps based on bam coverage files
##echo "computeMatrix reference-point --referencePoint TSS  -b 2000 -a 2000 -R TSS_heatmaps/$8 -S alignment_BAM/$1/$5\_shifted_FINAL.bw --missingDataAsZero --skipZeros -o TSS_heatmaps/$1/$5\_matrix.gz"
#->computeMatrix reference-point --referencePoint TSS  -b 2000 -a 2000 -R TSS_heatmaps/$8 -S alignment_BAM/$1/$5\_shifted_FINAL.bw --missingDataAsZero --skipZeros -o TSS_heatmaps/$1/$5\_matrix.gz

#echo "plotHeatmap -m TSS_heatmaps/$1/$5\_matrix.gz -out TSS_heatmaps/$1/$5\.pdf --colorList=\'white,blue\' --plotFileFormat pdf"
#->plotHeatmap -m TSS_heatmaps/$1/$5\_matrix.gz -out TSS_heatmaps/$1/$5\.pdf --colorList='white,blue' --plotFileFormat pdf

# Compute PCR bottlenecking coeffcient (PBC) 
# Input bam file used here is sorted bam after duplicates marking, and we need to sort bam by names
#->samtools sort -@ 10 -n -O BAM -o alignment_BAM/$1/$5\_dupMarked_TMP.bam alignment_BAM/$1/$5\_dupMarked.bam

#echo "bedtools bamtobed -bedpe -i alignment_BAM/$1/$5\_dupMarked_TMP.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v $7 | sort | uniq -c | awk'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > QC/$1/$5\_pbc.txt"
#->bedtools bamtobed -bedpe -i alignment_BAM/$1/$5\_dupMarked_TMP.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v $7 | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > QC/$1/$5\_pbc.txt

#->rm alignment_BAM/$1/$5\_dupMarked_TMP.bam

# This is to generate the plot for fragment size distribution
##echo "samtools view alignment_BAM/$1/$5\_shifted_FINAL.bam | awk \'\$9>0\' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's\/^\[ \\t\]*\/\/' > QC/$1/$5\_fragment_length_count.txt"
#->samtools view alignment_BAM/$1/$5\_shifted_FINAL.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > QC/$1/$5\_fragment_length_count.txt

# Calculate Fraction of reads in peaks (FRiP)
#->TOTALREADS=$(samtools view -c alignment_BAM/$1/$5\_shifted_FINAL.bam)
#->PEAKREADS=$(bedtools sort -i peaks/$1/$5\_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a alignment_BAM/$1/$5\_shifted_FINAL.bam  -b stdin -ubam | samtools view -c) 
#->echo $TOTALREADS $PEAKREADS $(awk "BEGIN {print $TOTALREADS/$PEAKREADS}") > QC/$1/$5\_FRiP.txt


