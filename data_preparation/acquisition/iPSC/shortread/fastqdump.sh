#!/usr/bin/env sh

# Author: Sean Maden
#
# Use fastq-dump to download a short-read RNA-seq sample's fastq's from SRA.
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

# paths
srrid=SRR6026510
outdir=./eternity/data/RI_benchmarking_hx1/

# download with fastq-dump
sudo ~/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump --accession $srrid --outdir $outdir --gzip