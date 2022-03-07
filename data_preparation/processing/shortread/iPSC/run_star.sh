#!/usr/bin/env sh

# Author: Sean Maden
#
# Run STAR and get STAR logs.

# star index path
srrid=SRR6026510
fastqdir=RI_benchmarking_fastqs
annodir=RI_benchmarking_resources/gencode_v35_annotation_files
STAR_INDEX_DIR=$annodir/STAR_index_gencode_v35
STAR_OUTPUT_DIR=RI_benchmarking_starlogs
READ1=$fastqdir/$srrid'_1.fastq.gz'
READ2=$fastqdir/$srrid'_2.fastq.gz'
OUTPREFIX=RI_benchmarking_starlogs/star_$srrid

sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode alignReads \
    --outSAMstrandField intronMotif --outFileNamePrefix  $OUTPREFIX \
    --genomeDir $STAR_INDEX_DIR --readFilesCommand zcat \
    --readFilesIn $READ1 $READ2

# count anno/non-anno junctions
# number of uniquely mapped reads: 53916311
# percent of uniquely mapped reads: 59.03%
cd RI_benchmarking_starlogs
awk -F '\t' '{print $6}' star_SRR6026510SJ.out.tab | sort | uniq -c | sort -nr
# 202939 1
# 105593 0
# 105593/(105593+202939) = 0.3422433 non-anno junctions

# compare with SRR2911306, hx1
# number of uniquely mapped reads: 21618226
# percent uniquely mapped reads: 88.37% 
awk -F '\t' '{print $6}' starSJ.out.tab | sort | uniq -c | sort -nr
# 162992 1
#  52885 0
# 52885/(52885+162992) = 0.2449775