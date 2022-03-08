#!/usr/bin/bash

# Authors: Mary Wood, Sean Maden
#
# UNSORTED_BAM = output unsorted bam file
# STAR_OUTPUT_DIR = output directory from STAR
# SORTED_BAM = output sorted bam file
# STAR_OUTPUT_DIR = output directory from STAR

# make unsorted bam
samtools view -b -o UNSORTED_BAM STAR_OUTPUT_DIR/Aligned.out.sam

# make sorted bam
samtools sort -O BAM -o SORTED_BAM STAR_OUTPUT_DIR/Aligned.out.sam

            