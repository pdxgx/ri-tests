#!/usr/bin/env

# Author: Sean Maden
# 
# Run IRFinder-S from a Docker image. Follow instructions in the GH repo
# to download the software and save its annotation files to the 
# `annodpath` variable. With the docker image we use -r and -d flags
# to specify the local paths to the intron annotations references and the
# output directory where results files will be saved.
# 
# Docker setup:
# 1. get the IRFinder-S image
# > sudo docker pull cloxd/irfinder:2.0
#
# 2. run the IRFinder-S image
# > sudo docker run cloxd/irfinder:2.0 --help
#

# manage paths
srrid=SRR2911306 # set run id
bamfpath=RI_benchmarking_BAMs/$srrid'.sorted.bam'
outdpath=RI_benchmarking_irfinders/
outrunpath=$outdpath$srrid/
annopath=RI_benchmarking_resources/gencode_v35_annotation_files/IRFinder_annotation

# navigate to main dir
cd /eternity/data/
# view docstrings
sudo docker run cloxd/irfinder:2.0 --help
# run irfinders using local annotation files
sudo docker run -w $PWD -v $PWD:$PWD cloxd/irfinder:2.0 -m BAM -r $refpath -d $outrunpath $bamfpath