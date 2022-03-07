#!/usr/bin/env R

# Author: Sean Maden
# 
# Filter and save long read BAM file

library(GenomicRanges)

srrid <- "SRR6026510"
run.handle <- "ipsc"

# long read granges
lr.gr.fname <- paste0("granges_longread_",srrid,".rda")
intv.gr <- get(load(lr.gr.fname))
intv.gr <- reduce(intv.gr)

library(Rsamtools)
bam.fname <- "SRP098984_SAMN07611993.merged.aligned.sorted.bam"
bam <- scanBam(bam.fname)