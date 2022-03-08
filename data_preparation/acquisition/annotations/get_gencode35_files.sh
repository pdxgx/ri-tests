#!/usr/bin/env bash

# Author Sean Maden
#
# Download human genome annotations in GTF and GFF3 formats, for genome
# build GENCODE v35, GRCh38/hg38. Also download the FASTA primary assembly 
# sequence.
#

urlbase=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/
fnv=(gencode.v35.primary_assembly.annotation.gtf.gz gencode.v35.annotation.gff3.gz GRCh38.primary_assembly.genome.fa.gz)

for fn in $fnv
	do fileurl=$urlbase$fn
		sudo wget $fileurl; sudo gunzip $fn
	done