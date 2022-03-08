#!/usr/bin/env sh

# Author: Sean Maden, Mary Wood
# 
# Prep a BAM file for retained introns analysis. Steps follow Mary's original 
# script. Code intended to be run in a remote server env where "sudo" access
# is usually required.

#----------------
# conda env setup
#----------------
conda create --name ri_bamprep
conda activate ri_bamprep
# check available versions
conda search -f star
conda search star
# install dependencies
conda install -c bioconda star=2.7.6a
conda install -c bioconda bowtie2=2.3.4.3
conda install -c bioconda samtools=1.3.1

#----------------
# download fastqs
#----------------
sudo /opt/anaconda3/envs/ri_bamprep/bin/fastq-dump --accession SRR6026510 \
    --split-3 --outdir . --gzip

#-----------------
# build star index
#-----------------
# build star index
sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode genomeGenerate \
    --genomeDir STAR_INDEX_DIR --genomeFastaFiles GENCODE_FASTA --sjdbGTFfile 

#-----------
# align BAMs
#-----------
# star index path
STAR_INDEX_DIR=RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/
STAR_OUTPUT_DIR=RI_benchmarking_hx1
READ1=RI_benchmarking_fastqs/SRR6026510.fastq
# align fastqs
sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode alignReads \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix RI_benchmarking_hx1/star \
    --genomeDir RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/ \
    --readFilesCommand zcat --readFilesIn RI_benchmarking_fastqs/SRR6026510_1.fastq.gz \
    RI_benchmarking_fastqs/SRR6026510_2.fastq.gz
# make unsorted bam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools view -b -o \
    RI_benchmarking_hx1/SRR6026510.unsorted.bam \
    RI_benchmarking_hx1/starAligned.out.sam
# make sorted bam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools sort -O BAM -o \
    RI_benchmarking_hx1/SRR6026510.sorted.bam \
    RI_benchmarking_hx1/SRR6026510.unsorted.bam