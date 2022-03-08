#!/usr/bin/env sh

# Author: Sean Maden
# 
# Prep a BAM file for retained introns analysis. Steps follow Mary's original 
# script. Code intended to be run in a remote server env where "sudo" access
# is usually required.

conda create --name ri_bamprep

conda activate ri_bamprep

# check available versions
conda search -f star
conda search star

# install dependencies
conda install -c bioconda star=2.7.6a
conda install -c bioconda bowtie2=2.3.4.3
conda install -c bioconda samtools=1.3.1

# get annotations

# DOWNLOAD FASTQ
# download
#sudo ~/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump --accession SRR2911306 \
# --outdir ./eternity/data/RI_benchmarking_hx1/ --split-files --gzip
sudo /opt/anaconda3/envs/ri_bamprep/bin/fastq-dump --accession SRR2911306 \
    --split-3 --outdir . --gzip
# returns
# Read 24463210 spots for SRR2911306
# Written 24463210 spots for SRR2911306
# expand the gz file
# sudo gzip -d SRR2911306.fastq.gz

# PROCESS BAM FILES
#gtfpath=RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gtf
#fapath=RI_benchmarking_resources/gencode_v35_annotation_files/GRCh38.primary_assembly.genome.fa
# build star index
# STAR --runMode genomeGenerate --genomeDir STAR_INDEX_DIR --genomeFastaFiles GENCODE_FASTA --sjdbGTFfile 

# star index path
STAR_INDEX_DIR=RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/
STAR_OUTPUT_DIR=RI_benchmarking_hx1
READ1=RI_benchmarking_hx1/SRR2911306.fastq

# align fastq
#STAR --runMode alignReads --outSAMstrandField intronMotif --outFileNamePrefix STAR_OUTPUT_DIR \
#                --genomeDir STAR_INDEX_DIR --readFilesCommand zcat --readFilesIn READ1 READ2

sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode alignReads \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix RI_benchmarking_hx1/star \
    --genomeDir RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/ \
    --readFilesCommand zcat --readFilesIn RI_benchmarking_hx1/SRR2911306_1.fastq.gz \
    RI_benchmarking_hx1/SRR2911306_2.fastq.gz
# returns
# Nov 10 18:37:20 ..... started STAR run
# Nov 10 18:37:20 ..... loading genome
# Nov 10 18:37:52 ..... started mapping
# Nov 10 19:12:04 ..... finished mapping
# Nov 10 19:12:08 ..... finished successfully

# make unsorted bam
# samtools view -b -o UNSORTED_BAM STAR_OUTPUT_DIR/Aligned.out.sam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools view -b -o \
    RI_benchmarking_hx1/SRR2911306.unsorted.bam \
    RI_benchmarking_hx1/starAligned.out.sam
# make sorted bam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools sort -O BAM -o \
    RI_benchmarking_hx1/SRR2911306.sorted.bam \
    RI_benchmarking_hx1/SRR2911306.unsorted.bam

# PREP THE KMA RANGES
seqpath=seq.fa
gtfpath=trans.gtf
outdir=RI_benchmarking_hx1
PRE=/path/to/kma/preprocess
python $PRE/generate_introns.py --genome seq.fa --gtf trans.gtf --extend N --out out_dir






