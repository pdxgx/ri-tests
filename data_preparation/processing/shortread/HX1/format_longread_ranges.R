#!/usr/bin/env R

# Author: Sean Maden
# Format the table of long read ranges.

# get the long read ranges
tsv.fname <- "RI_txs_to_read_ids_final_09-21-2021_16.53.14.tsv"
tsv <- fread(tsv.fname, sep = "\t", header = T, data.table = F)

# long read intron granges
intv <- unlist(lapply(unique(tsv$transcript), function(tid){
  message(tid); tsvf <- tsv[which(tsv$transcript == tid)[1],]
  intv <- unique(unlist(strsplit(tsvf$tx_introns, ";")))
  intv <- paste0(tsvf$chrom[1], ":", intv); return(intv)}))
df.intv <- data.frame(chr = gsub(":.*", "", intv),
                      start = gsub(".*:|-.*", "", intv),
                      end = gsub(".*-", "", intv), 
                      stringsAsFactors = F)
intv.gr <- GenomicRanges::makeGRangesFromDataFrame(df.intv)

# save
save(intv.gr, file = "granges_longread_SRR6026510.rda")