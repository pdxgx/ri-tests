#!/usr/bin/env bash

# Author: Sean Maden
# 
# Run the iREAD software to quantify intron expression from 
# STAR-aligned BAMs.
# 
# NOTES:
# 
# * get gff annotation file:
#       -- https://www.gencodegenes.org/human/
#       -- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz
# ENV:
# 
# activate the conda env with python 2###
# > conda activate py2718_iread

#--------------------
# manage vars & paths
#--------------------
srrid=SRR6026510
bamdpath='/eternity/data/RI_benchmarking_BAMs'
bedpath='/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.introns.bed' 
ireadpath=iread.py #/home/metamaden/iread/iread.py
outputpath=$srrid'_iread'
sudo mkdir $outputpath

#----------
# run iread
#----------
cd iread
sudo chmod -R 777 .
# iterate on chr chunks, implicit sudo access
for chri in $chrv
    do
        echo 'Beginning chr '$chri
        bamfpath=$bamdpath/$srrid'.sorted.'$chri'.nochr.bam'
        numreads=$(samtools view -c -F 260 $bamfpath)
        python iread.py $bamfpath $bedpath -o $outputpath -t $numreads
        echo 'Finished chr '$chri
    done