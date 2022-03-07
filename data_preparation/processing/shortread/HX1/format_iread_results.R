#!/usr/bin/env R

# Author: Sean Maden
# Format/clean up the iREAD results, removing duplicated regions 
# with null results. Then get the data formatted as a GRanges object.

library(dplyr)
library(GenomicRanges)

srrid <- "SRR2911306"; run.handle <- "hx1"

#----------
# load data
#----------
ir.fname <- paste0(srrid, "_allchr_iread.txt")
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
length(dup.id) 
# 242035

# get new df for duplicated ids
which.notdup <- ir$intronid %in% dup.id & !duplicated(ir$intronid)
df.notdup <- ir[which.notdup,]
dim(df.notdup) 
# [1] 242035     11
# retain max observed expr by dup id
dfply <- ir[ir$intronid %in% dup.id,]; dim(dfply) 
# 968140     11
eval.str <- paste0("dfply <- dfply %>% group_by(intronid) ",
                   "%>% summarise(",cname," = max(",cname,", na.rm = T))")
eval(parse(text = eval.str)); dim(dfply) 
# 242035      2
dfply <- dfply[order(match(dfply$intronid, df.notdup$intronid)),]
identical(dfply$intronid, df.notdup$intronid) 
# TRUE
df.notdup[,cname] <- dfply[,2]

# bind results
irf <- rbind(ir[!ir$intronid %in% dup.id,], df.notdup)
dim(ir) 
# 968140     11
dim(irf) 
# 242035     11

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
