#!/usr/bin/env sh

# Author: Sean Maden
# 
# Download the short-read RNA-seq FASTQ files from SRA.
srrid=SRR6026510
sudo /opt/anaconda3/envs/ri_bamprep/bin/fastq-dump --accession $srrid \
	--split-3 --outdir . --gzip