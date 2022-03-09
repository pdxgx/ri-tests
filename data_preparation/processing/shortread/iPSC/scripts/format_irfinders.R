#!/usr/bin/env R

# Author: Sean Maden
#
# Format results output from IRFinder-S prior to LWM harmonization.
# 

# set run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

#----------
# load data
#----------
irfs.fname <- "IRFinder-IR-nondir.txt"
irfs <- data.table::fread(irfs.fname, header = T, data.table = F)

#--------
# granges
#--------
# make granges
colnames(irfs)[1:3] <- c("chr", "start", "end")
irfs.gr <- GenomicRanges::makeGRangesFromDataFrame(irfs, keep.extra.columns = T)

# save the granges
irfs.gr.fname <- paste0("granges_irfinders_",
                       srrid,"-",run.handle,".rda")
save(irfs.gr, file = irfs.gr.fname)

# save the df
irfs.df <- as.data.frame(irfs.gr, stringsAsFactors = F)
irfs.df.fname <- paste0("df-granges_irfinders_",
                       srrid,"-",run.handle,".csv")
write.csv(irfs.df, file = irfs.df.fname)