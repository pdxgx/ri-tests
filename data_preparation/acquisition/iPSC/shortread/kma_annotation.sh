#!/usr/bin/env sh

# Authors: Sean Maden, Mary Wood
#
# Make intron annotations for KMA. Make an intron annotations
# sequence FASTQ using KMA's Python script and a GENCODE v35 GTF. 
# Run this script from a clone of the main KMA repo.
#
# This script saves 3 new annotation files: 
# > introns.bed
# > introns.fa
# > and intron_to_transcripts.txt

# get paths
pyscriptfpath=/usr/lib64/R/library/kma/pre-process/generate_introns.py
annopath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/
fastafname=GRCh38.primary_assembly.genome.fa
gtffname=gencode.v35.annotation.gtf
outdname=kma_annotations
fapath=$annopath$fastafname
gtfpath=$annopath$gtffname
outdirpath=$annopath$outdname

# generate the annotation files
sudo python2 $preprocesspath/$prescript --genome $fapath --gtf $gtfpath --extend 25 \
	--out $outdirpath

# combine the introns and transcripts `fa` files
sudo bash -c 'cat ./kma_reference/gencode.v35.transcripts.fa ./sm/introns.fa > trans_and_introns_sm.fa'
