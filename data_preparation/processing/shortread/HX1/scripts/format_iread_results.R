#!/usr/bin/env R

# Author: Sean Maden
# 
# Format/clean up the iREAD results prior to LWM harmonization. First, remove 
# duplicated regions with null results. Then get the data formatted as a 
# GRanges object.
#

library(dplyr)
library(GenomicRanges)

# set run identifier
srrid <- "SRR2911306"
run.handle <- "hx1"

#----------
# load data
#----------
ir.fname <- paste0(srrid,"_allchr_iread.txt")
ir <- read.table(ir.fname, sep = " ", header = T)

#-------------------------
# filter duplicate regions
#-------------------------
cname <- "fpkm"

# get the feature labels/id
list.intron <- strsplit(ir$intron_id, "-")
ir$chr <- paste0("chr", unlist(lapply(list.intron, function(i){i[[1]]})))
ir$start <- unlist(lapply(list.intron, function(i){i[[2]]}))
ir$end <- unlist(lapply(list.intron, function(i){i[[3]]}))
ir$intronid <- paste0(ir$chr, ":", ir$start, "-", ir$end)

# get duplicated features
dup <- which(duplicated(ir$intronid))
dup.id <- unique(ir$intronid[dup])

# get new df for duplicated ids
which.notdup <- ir$intronid %in% dup.id & !duplicated(ir$intronid)
df.notdup <- ir[which.notdup,]
# retain max observed expr by dup id
dfply <- ir[ir$intronid %in% dup.id,]
eval.str <- paste0("dfply <- dfply %>% group_by(intronid) ",
                   "%>% summarise(",cname," = max(",cname,", na.rm = T))")
eval(parse(text = eval.str))
dfply <- dfply[order(match(dfply$intronid, df.notdup$intronid)),]
df.notdup[,cname] <- dfply[,2]
# bind results
irf <- rbind(ir[!ir$intronid %in% dup.id,], df.notdup)

#--------
# granges
#--------
irf.gr <- makeGRangesFromDataFrame(irf, keep.extra.columns = T)

# save the granges
irf.gr.fname <- paste0("granges_iread_",
                       srrid,"-",run.handle,".rda")
save(irf.gr, file = irf.gr.fname)

# save the df
irf.df <- as.data.frame(irf.gr, stringsAsFactors = F)
irf.df.fname <- paste0("df-granges_iread_",
                       srrid,"-",run.handle,".csv")
write.csv(irf.df, file = irf.df.fname)