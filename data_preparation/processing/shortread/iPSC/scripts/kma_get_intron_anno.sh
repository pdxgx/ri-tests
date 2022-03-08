#!/usr/bin/env sh

# Author: Sean Maden
#
# Get intron ranges from GTF file using KMA preprocess script.
#

# download the grch38 fasta
faurl=http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget $fwurl

# make symlink to file
fapath=/Users/maden/Desktop/misc/retained_intron_RI/Homo_sapiens.GRCh38.dna.primary_assembly.fa
sypath=/Users/maden/Desktop/misc/retained_intron_RI/SRR2911306_hx1/gencode38/GRCh38.dna.primary_assembly.fa
ln -s $fapath $sypath

# install kma R package, get path to scripts
conda install -c bioconda kma
kmapath=/Library/Frameworks/R.framework/Versions/4.1/Resources/library/kma/pre-process
scriptpath=$kmapath/generate_introns.py

# get paths to fasta, gtf
gtffname=gencode.v38.nochr.annotation.gtf
gtfpath=/Users/maden/Desktop/misc/retained_intron_RI/SRR2911306_hx1/gencode38/$gtffname

# run kma preprocess to get introns gtf
nvar=25 # amount to extend
python2 $scriptpath --genome $sypath --gtf $gtfpath --extend $nvar --out .
