#!/usr/bin/env R

# Author: Sean Maden
# 
# Format the MAJIQ results and store as GenomicRanges object.
#

library(GenomicRanges)

#----------
# load data
#----------
# set the run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

# load the majiq results
rm.fname <- paste0(run.handle, "_RI.MATS.JC.TXT")
rm.dir <- paste0(srrid, "_", run.handle)
rm.fpath <- file.path(rm.dir, rm.fname)
rm <- read.table(rm.fpath, sep = "\t", header = T)

# remove na inclusion levels
dim(rm) # [1] 7317   23
rm <- rm[!is.na(rm$IncLevel1),]
dim(rm) # [1] 6093   23

#--------------------
# get unique features
#--------------------
coordv <- paste0(rm$chr, ":", rm$riExonStart_0base, "-", rm$riExonEnd)
dft <- as.data.frame(table(coordv))
max(dft[,2]) # [1] 8

# get unique introns
rm$coord.str <- coordv
coordv.unique <- dft[dft[,2]==1,1]
dfir <- rm[rm$coord.str %in% coordv.unique,]
dfir <- data.frame(coord = dfir$coord.str, rmats.inclvl = dfir$IncLevel1)
dim(dfir) # [1] 4184    2

# get collapsed inclusions for repeated features
coordv.dup <- as.character(dft[dft[,2]>1,1])
dfir <- rbind(dfir, do.call(rbind, lapply(coordv.dup, function(ir.coordi){
  dfii <- rm[rm$coord.str==as.character(ir.coordi),]
  medi <- median(as.numeric(dfii$IncLevel1), na.rm = T)
  data.frame(coord = ir.coordi, rmats.inclvl = medi)
})))
dim(dfir)
# [1] 5033    2

#--------------------
# make and save grset
#--------------------
# separate coords
dfir$chr <- gsub(":.*", "", dfir$coord)
dfir$start <- gsub("-.*|.*:", "", dfir$coord)
dfir$stop <- gsub(".*-", "", dfir$coord)
rm.gr <- makeGRangesFromDataFrame(dfir[,c(2:5)], keep.extra.columns = T)

# save granges
rm.gr.fname <- paste0("granges_majiq_", srrid,"-",run.handle,".rda")
rm.gr.fpath <- file.path(rm.dir, rm.gr.fname)
save(rm.gr, file = rm.gr.fpath)