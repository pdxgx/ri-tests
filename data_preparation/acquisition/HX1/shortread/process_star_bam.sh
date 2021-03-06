#!/usr/bin/env sh

# Author: Sean Maden, Mary Wood
# 
# Prep a STAR-aligned, sorted BAM file for retained introns analysis. 
# Code intended to be run in a remote server env where "sudo" access is 
# usually required. 
# 
# Code steps include: 
#   (1.) env setup using conda; 
#   (2.) fastq download with sra-toolkit/fastq-dump; 
#   (3.) STAR index creation; 
#   (4.) BAM alignment with STAR index; 
#   (5.) BAM sorting with samtools.
#
# To setup the sra-toolkit, use either:
# 
# 1. conda env
# > conda activate sra_process
# > conda install -c bioconda/label/cf201901 sra-tools
# 
# 2. download locally
# > wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-centos_linux64.tar.gz
# > tar xzvf sratoolkit.2.9.6-centos_linux64.tar.gz
# > export PATH=$PATH'/home/metamaden/sratoolkit.2.9.6-centos_linux64/bin'

#------------------------------
# (0.) get file names and paths
#------------------------------
GENCODE_FASTA= # path to gencode fasta
srrid=SRR2911306
outdpath=RI_benchmarking_hx1/
STAR_INDEX_DIR=RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/
STAR_OUTPUT_DIR=$outdpath
outfileprefix=$outdpath'star'
starindexdpath=RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35/
fastqpath1=RI_benchmarking_fastqs/$srrid'_1.fastq.gz'
fastqpath2=RI_benchmarking_fastqs/$srrid'_2.fastq.gz'
# sam and bam paths
samfpath=$outdpath'starAligned.out.sam'
unsortedbamfpath=$outdpath$srrid'.unsorted.bam'
sortedbamfpath=$outdpath$srrid'.sorted.bam'

#---------------------
# (1.) conda env setup
#---------------------
conda create --name ri_bamprep
conda activate ri_bamprep
# check available versions
conda search -f star
conda search star
# install dependencies
conda install -c bioconda star=2.7.6a
conda install -c bioconda bowtie2=2.3.4.3
conda install -c bioconda samtools=1.3.1

#---------------------
# (2.) download fastqs
#---------------------
sudo /opt/anaconda3/envs/ri_bamprep/bin/fastq-dump --accession $srrid \
    --split-3 --outdir . --gzip

#----------------------
# (3.) build star index
#----------------------
# build star index
sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode genomeGenerate \
    --genomeDir STAR_INDEX_DIR --genomeFastaFiles GENCODE_FASTA --sjdbGTFfile 

#------------------
# (4.) align fastqs
#------------------
sudo /opt/anaconda3/envs/ri_bamprep/bin/STAR --runMode alignReads \
    --outSAMstrandField intronMotif --outFileNamePrefix $outfileprefix \
    --genomeDir $starindexdpath --readFilesCommand zcat \
    --readFilesIn $fastqpath1 $fastqpath2

#---------------
# (5.) sort bams
#---------------
# make unsorted bam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools view -b -o \
    $unsortedbamfpath $samfpath
# make sorted bam
sudo /opt/anaconda3/envs/ri_bamprep/bin/samtools sort -O BAM -o \
    $sortedbamfpath $unsortedbamfpath