#!/usr/bin/env R

# Author: Sean Maden
#
# Format results output from IRFinder-S

srrid <- "SRR2911306"; run.handle <- "hx1"

#----------
# load data
#----------
irfs.fname <- "IRFinder-IR-nondir.txt"
irfs <- data.table::fread(irfs.fname, header = T, data.table = F)

table(irfs$Warnings)
#
#     -              LowCover           LowSplicing          MinorIsoform NonUniformIntronCover 
# 59604                161017                  2227                  6317                 20814

table(irfs$Warnings)/nrow(irfs)
#           -              LowCover           LowSplicing          MinorIsoform NonUniformIntronCover 
# 0.238436029           0.644122106           0.008908748           0.025270123           0.083262994

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