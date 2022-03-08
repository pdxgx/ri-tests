#!/usr/bin/env R

# Author: Sean Maden
#
# Do samtools mpileup, getting summaries by gene.

library(rtracklayer)

#-----------
# load data
#-----------
bed.fname <- "gencode.v38.annotation.bed"
bed <- read.table(bed.fname)


# identify genes from gencode gtf
bedname_read=gencode.v38.annotation.bed
bedname_write=gencode_genes.v38.annotation.bed
# get gene bed file
awk '$8 == "gene"' $bedname_read > $bedname_write

# define vars
bamfname='SRR6026510.sorted.bam'; mpiletmp='mpiletemp'
bedtmp='bedtmp.bed'

# test 1st line
line=$(head -n 1 $bedname_write)
head -n 1 $bedname_write > $bedtmp

line=$p
# get gene point info
chri=$(cut -f1 $bedtmp) 
geneidi=$(cut -f4 $bedtmp)
# get gene length
coordstart=$(cut -f2 $bedtmp); coordend=$(cut -f3 $bedtmp)
coorddiff=$((coordstart-coordend)); genelength=${coorddiff#-}

# run samtools, get coverage/bp for the gene
samtools mpileup -r $chri -l $bedtmp $bamfname > $mpiletmp