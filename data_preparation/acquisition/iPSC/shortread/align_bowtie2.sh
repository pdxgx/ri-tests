#!/usr/bin/env sh

# Authors: Sean Maden, Mary Wood
#
# Align FASTQs using bowtie2 and KMA intron annotation file. Run this
# script after running `make_kma_annotation.sh` to generate the inton
# annotation file.

srrid=SRR6026510

# prepare bt2 file and new sample bam files
sudo bowtie2-build --offrate 1 trans_and_introns.fa trans_and_introns_new

# get file paths
fqfpath1='/eternity/data/RI_benchmarking_fastqs/'$srrid'_1.fastq'
fqfpath2='/eternity/data/RI_benchmarking_fastqs/'$srrid'_2.fastq'
tipath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_new
newbampath=/home/metamaden/newbam.bam

# run alignment 
sudo bowtie2 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x $tipath -1 \
$fqfpath1 -2 $fqfpath2 | samtools view -Sb - > $newbampath