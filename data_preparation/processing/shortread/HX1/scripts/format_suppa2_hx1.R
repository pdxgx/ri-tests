#!/usr/bin/env R

# Author: Sean Maden
#
# Format SUPPA2 outputs for harmonization.

library(GenomicRanges)

#----------
# load data
#----------
# set the run identifiers
srrid <- "SRR2911306"
run.handle <- "hx1"

# load the suppa2 results
sp.fname <- "hx1_ri-tpm-suppa2.txt.psi"
sp.dir <- paste0(srrid, "_", run.handle)
sp.fpath <- file.path(sp.dir, sp.fname)
sp <- read.table(sp.fpath, sep = "\t", header = T)
sp <- data.frame(transcript_id=rownames(sp), psi=as.numeric(sp[,1]))

# remove na inclusion levels
dim(sp) # [1] 7540    2
sp <- sp[!is.na(sp$psi),]
dim(sp) # [1] 5875    2

# check for unique features
length(unique(sp$transcript_id)) # 5875

#------------------------
# make and save new grset
#------------------------
# assign coordinates
coordv <- gsub(".*RI:", "", sp$transcript_id)
splitv <- strsplit(coordv, ":|-")
chrv <- sapply(splitv, function(ii){ii[1]})
startv <- sapply(splitv, function(ii){
  min(as.numeric(ii[2:5]))}) # max of coord
endv <- sapply(splitv, function(ii){
  max(as.numeric(ii[2:5]))}) # max of coord

# new grset
sp.df <- data.frame(chr = chrv, start = startv, end = endv,
                    psi = sp$psi)
sp.gr <- makeGRangesFromDataFrame(sp.df, keep.extra.columns = T)

# save granges
sp.gr.fname <- paste0("granges_suppa2_", srrid,"-",run.handle,".rda")
sp.gr.fpath <- file.path(sp.dir, sp.gr.fname)
save(sp.gr, file = sp.gr.fpath)