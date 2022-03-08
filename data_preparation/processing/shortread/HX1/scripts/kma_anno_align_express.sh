#!/bin/bash

# Authors: Sean Maden, Mary Wood
# 
# Steps to prep and run KMA. First makes `trans_and_introns.fa` file using
# Python script "generate_introns.py". Second, runs bowtie2 alignments using 
# KMA's custom intron annotations. Third, performs KMA on the aligned BAMs.
#

#-----------------------------
# 1. Make trans_and_introns.fa
#-----------------------------
# this populates the output dir with 3 files: 
# introns.bed, introns.fa, and intron_to_transcripts.txt
preprocesspath=/usr/lib64/R/library/kma/pre-process
prescript=generate_introns.py
annopath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files
fapath=$annopath/GRCh38.primary_assembly.genome.fa
gtfpath=$annopath/gencode.v35.annotation.gtf
outdir=$annopath/sm
sudo python2 $preprocesspath/$prescript --genome $fapath --gtf $gtfpath --extend 25 --out $outdir
# combine the introns and transcripts `fa` files
sudo bash -c 'cat ./kma_reference/gencode.v35.transcripts.fa ./sm/introns.fa > trans_and_introns_sm.fa'

#---------------------
# 2. align with bowtie
#---------------------
# prepare bt2 file and new sample bam files
sudo bowtie2-build --offrate 1 trans_and_introns_sm.fa trans_and_introns_sm
# get file paths
tipath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_sm
fqfpath1=/eternity/data/RI_benchmarking_hx1/SRR2911306_1.fastq
fqfpath2=/eternity/data/RI_benchmarking_hx1/SRR2911306_2.fastq
newbampath=/home/metamaden/newbam2.bam

# run alignment 
sudo bowtie2 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x $tipath -1 $fqfpath1 -2 $fqfpath2 | samtools view -Sb - > $newbampath

sudo bowtie2 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x $tipath -1 $fqfpath1 -2 $fqfpath2

# sudo gzip $fqfpath1; sudo rm $fqfpath2
# sudo rm $fqfpath1; sudo rm $fqfpath2

#---------------
# 3. run express
#---------------
# Quantify intron expression with eXpress. 
# This makes the sample `*.kma.bam` files

# add express to path
export PATH="$HOME/express-1.5.1-linux_x86_64:$PATH"

# get file paths
resultspath=/eternity/data/RI_benchmarking_results; outdir=express_output_sm
outpath=$resultspath/$outdir; ripath=/eternity/data/RI_benchmarking_BAMs
bamfname=SRR1660804_sm.kma.bam; bamfpath=$ripath/$bamfname
annodpath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files
fapth=$annodpath/trans_and_introns_sm.fa
sudo mkdir $outpath/$s
sudo ~/express-1.5.1-linux_x86_64/express -o $outpath/$s/ $fapth $bamfpath
samples=`ls /$resultspath/$olddir`
cd $ripath; samples=`ls -1 *$sjoutext  | sed -E 's/(.+)[.]SJ[.]out[.]tab/\1/'`
for s in $samples
    do echo $s
        sudo mkdir $outpath/$s
        sudo $expresscall -o $outpath/$s/ $fapth $bampath/$s$kmaext
    done

#fapth=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns.fa
#bampath=/eternity/data/RI_benchmarking_BAMs; tifname=trans_and_introns_sm.fa
# expresscall=express # /eternity/data/RI_benchmarking_resources/express-1.5.1-linux_x86_64/express
# fapth=$annodpath/$tifname; kmaext=.kma.bam; sjoutext=.SJ.out.tab
