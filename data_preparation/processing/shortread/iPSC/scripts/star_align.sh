#!/usr/bin/bash

# Authors: Mary Wood, Sean Maden
#
# STAR_OUTPUT_DIR = output directory
# STAR_INDEX_DIR = directory with STAR index
# READ1/READ2 = paired end FASTQ files (compressed with gzip)

STAR --runMode alignReads --outSAMstrandField intronMotif --outFileNamePrefix \
    STAR_OUTPUT_DIR --genomeDir STAR_INDEX_DIR --readFilesCommand zcat \
    --readFilesIn READ1 READ2

