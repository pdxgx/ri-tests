#!/usr/bin/env bash

# Author: Sean Maden
# Script to prepare KMA run

# Notes:
# * genome seq fasta: download latest hg38 version from: 
#        -- https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39
#        -- https://www.ncbi.nlm.nih.gov/genome/guide/human/
#        -- http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

conda activate r351_kma

#-------------------
# get example fastqs
#-------------------
fastq-dump -I --split-files SRR5009474
fastq-dump -I --split-files SRR1660804


#---------------
# get main paths
#---------------

kmadir=/Users/maden/miniconda3/envs/r351_kma/lib/R/library/kma
scriptdir=pre-process
annodir=./genome_anno

#----------
# file prep
#----------
# convert `chr##` to `##` in gtf file
gtf=$annodir/gencode.v38.annotation.gtf
gtfnochr=$annodir/gencode.v38.nochr.annotation.gtf
sed 's/chr//g' $gtf > $gtfnochr

#--------------
# preprocessing 
#--------------
# notes: generates intron coordinates from seq fasta and gtf anno

PRE=$kmadir/$scriptdir
seq=$annodir/hg38.fa # note: chr indices need to match those in gtf, e.g. "chr1"
gtf=$annodir/gencode.v38.annotation.gtf # note: chr indices need to match those in *.fa file, e.g. "chr1"
outdir=kma_output
nvar=25

python $PRE/generate_introns.py --genome $seq --gtf $gtf --extend $nvar --out $outdir

# returns: 
# INFO: Reading in GTF: ./genome_anno/gencode.v38.annotation.gtf
# INFO: Grouping transcripts by gene
# INFO: Writing intron BED file: kma_output/introns.bed
# INFO: Computing intron-to-transcript compatability
# INFO: Opening FASTA: ./genome_anno/hg38.fa
# INFO: Note: will take a while the first time it is opened.
# INFO: On intron: 10000
# INFO: On intron: 20000
# INFO: On intron: 30000
# INFO: On intron: 40000
# INFO: On intron: 50000
# INFO: On intron: 60000
# INFO: On intron: 70000
# INFO: On intron: 80000
# INFO: On intron: 90000
# INFO: On intron: 100000
# INFO: On intron: 110000
# INFO: On intron: 120000
# INFO: On intron: 130000
# INFO: On intron: 140000
# INFO: On intron: 150000
# INFO: On intron: 160000
# INFO: On intron: 170000
# INFO: On intron: 180000
# INFO: On intron: 190000
# INFO: On intron: 200000
# INFO: Writing intron sequences out to kma_output/introns.fa

# combine transcripts and introns fasta files
# > cat trans.fa introns.fa > trans_and_introns.fa
cat $seq $outdir/introns.fa > $outdir/trans_and_introns.fa

# make the bowtie2 index
bowtie2-build --offrate 1 $outdir/trans_and_introns.fa $outdir/trans_and_introns

# returns:
# Settings:
#   Output files: "kma_output/trans_and_introns.*.bt2"
#   Line rate: 6 (line is 64 bytes)
#   Lines per side: 1 (side is 64 bytes)
#   Offset rate: 1 (one in 2)
#   FTable chars: 10
#   Strings: unpacked
#   Max bucket size: default
#   Max bucket size, sqrt multiplier: default
#   Max bucket size, len divisor: 4
#   Difference-cover sample period: 1024
#   Endianness: little
#   Actual local endianness: little
#   Sanity checking: disabled
#   Assertions: disabled
#   Random seed: 0
#   Sizeofs: void*:8, int:4, long:8, size_t:8
# Input files DNA, FASTA:
#   kma_output/trans_and_introns.fa
# Building a SMALL index
# Reading reference sizes
#   Time reading reference sizes: 00:00:53
# Calculating joined length
# Writing header
# Reserving space for joined string
# Joining reference sequences
#   Time to join reference sequences: 00:00:36
# bmax according to bmaxDivN setting: 1010782202
# Using parameters --bmax 758086652 --dcv 1024
#   Doing ahead-of-time memory usage test
#   Passed!  Constructing with these parameters: --bmax 758086652 --dcv 1024
# Constructing suffix-array element generator
# Building DifferenceCoverSample
#   Building sPrime
#   Building sPrimeOrder
#   V-Sorting samples
#   V-Sorting samples time: 00:03:04
#   Allocating rank array
#   Ranking v-sort output
#  Ranking v-sort output time: 00:00:42
#  Invoking Larsson-Sadakane on ranks
#  Invoking Larsson-Sadakane on ranks time: 00:00:58
#  Sanity-checking and returning
# Building samples
# Reserving space for 12 sample suffixes
# Generating random suffixes
# QSorting 12 sample offsets, eliminating duplicates
# QSorting sample offsets, eliminating duplicates time: 00:00:00
# Multikey QSorting 12 samples
#   (Using difference cover)
#   Multikey QSorting samples time: 00:00:00
# Calculating bucket sizes
# Splitting and merging
#   Splitting and merging time: 00:00:00
# Avg bucket size: 4.04313e+09 (target: 758086651)
# Converting suffix-array elements to index image
# Allocating ftab, absorbFtab
# Entering Ebwt loop
# Getting block 1 of 1
#   No samples; assembling all-inclusive block
#   Sorting block of length 4043128809 for bucket 1
#   (Using difference cover)

# align reads for a single paired read sample
fq1=samples/fastq/SRR1660804_1.fastq
fq2=samples/fastq/SRR1660804_2.fastq

bowtie2 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -X trans_an
 d_introns
-1 left.fastq -2 right.fastq | samtools view -Sb - > hits.bam




