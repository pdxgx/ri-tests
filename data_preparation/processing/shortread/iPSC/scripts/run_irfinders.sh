#!/usr/bin/env

# Author: Sean Maden

# Run IRFinder-S
# 
# Docker setup
# a. get the image
# > sudo docker pull cloxd/irfinder:2.0
# b. run the image
# > sudo docker run cloxd/irfinder:2.0 --help

# view docstrings
sudo docker run cloxd/irfinder:2.0 --help

# run using local refs
cd /eternity/data/
outdpath=RI_benchmarking_irfinders
refpath=RI_benchmarking_resources/gencode_v35_annotation_files/IRFinder_annotation

# run for SRR6026510 ipsc
bamfpath=RI_benchmarking_BAMs/SRR6026510.sorted.bam
outrunpath=$outdpath/SRR6026510_ipsc/
sudo docker run -w $PWD -v $PWD:$PWD cloxd/irfinder:2.0 -m BAM -r $refpath -d $outrunpath $bamfpath
