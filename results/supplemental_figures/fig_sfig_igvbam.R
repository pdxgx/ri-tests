#!/usr/bin/env R

# Author: Sean Maden
#
# Get the filtered reads for the iPSC LR BAM.

library(Rsamtools)

# set the tsv fpath
tsv.fpath <- file.path("RI_txs_to_read_ids_final_01-06-2022_17.05.22_iPSC.tsv")
tsv <- data.table::fread(tsv.fpath)

# filter bam on qname (read ID column)
read.idv <- tsv$orig_read_id
newbam.fname <- "lr-ipsc_filtered-reads.bam"
filter <- FilterRules(list(samplefilt=function(x){x$qname %in% read.idv}))
filterBam(bam.fname, newbam.fname, filter=filter)

