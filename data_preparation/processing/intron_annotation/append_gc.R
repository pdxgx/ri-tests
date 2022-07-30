#!/usr/bin/env R

# Author: Sean Maden
#
# Append GC content to intron annotations tables.

library(GenomicRanges)
library(Biostrings)
library(data.table)
# BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome"))
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)

# file file names vector
fnv <- c("nonzero_RI_data_summary_iPSC.tsv",
         "called_RI_data_summary_iPSC.tsv",
         "called_RI_data_summary_HX1.tsv",
         "nonzero_RI_data_summary_HX1.tsv")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# append and resave files with gc content
for(fni in fnv){
  message("beginning file ", fni, "...")
  # load data
  df <- data.table::fread(fni, sep = "\t", header = T, 
                          data.table = F)
  # get introns genomic ranges
  intv <- df$intron
  int.df <- data.frame(chr = gsub(":.*", "", intv),
                       start = gsub(".*:|-.*", "", intv),
                       end = gsub(".*-", "", intv))
  int.gr <- makeGRangesFromDataFrame(int.df)
  # get seqs 
  seqs <- BSgenome::getSeq(bsgenome, int.gr)
  # apply gc content
  lfv <- letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)
  # append data, write new table
  df$gc.fract <- as.numeric(lfv)
  new.fni <- paste0(gsub("\\.tsv$", "", fni), "_gc.tsv")
  data.table::fwrite(df, file = new.fni, sep = "\t", row.names = F)
  message("finished with file ", fni, ".")
}