#!/usr/bin/bash

# Authors: Mary, Sean

#STAR_INDEX_DIR = directory to hold STAR index
#GENCODE_FASTA = reference genome FASTA file from GENCODE
#GENCODE_GTF = GTF file from GENCODE

STAR --runMode genomeGenerate --genomeDir STAR_INDEX_DIR --genomeFastaFiles \
    GENCODE_FASTA --sjdbGTFfile GENCODE_GTF

