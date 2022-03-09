#!/usr/bin/env R

# Author: Sean Maden
# 
# Format the IntEREst results prior to LWM harmonization.
#

library(dplyr)
library(GenomicRanges)

# set run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

#----------
# load data
#----------
# load the results table
it.fname <- paste0(srrid, "_allchr_interest.tsv")
it <- read.table(it.fname, sep = "\t")
it <- it[it$int_ex == "intron",]
# set the colnames
colnames(it) <- c("chr", "start", "end", "strand",
                  "int_ex", "int_ex_num", "collapsed_transcripts_id",
                  "collapsed_transcripts", "collapsed_gene_id",
                  "IntRet_frequency", "IntRet_genewide_scaled")

#---------------------------------------
# filter max expr for duplicated introns
#---------------------------------------
cname <- "IntRet_frequency" # expr col of interest
it.df <- as.data.frame(it, stringsAsFactors = F)
it.df$intronid <- paste0(it.df$chr, ":", it.df$start, "-", it.df$end)

# get the duplicated ids
dup <- which(duplicated(it.df$intronid))
dup.id <- unique(it.df$intronid[dup])

#---------------------
# get granges and save
#---------------------
it <- it.df # it.dff
it.gr <- makeGRangesFromDataFrame(it, keep.extra.columns = T)

# save granges
it.gr.fname <- paste0("granges_interest_",
                      srrid,"-",run.handle,".rda")
save(it.gr, file = it.gr.fname)

# save df
it.gr.df <- as.data.frame(it.gr, stringsAsFactors = FALSE)
it.gr.df.fname <- paste0("df-granges_interest_",
                         srrid,"-",run.handle,".csv")
write.csv(it.gr.df, file = it.gr.df.fname)