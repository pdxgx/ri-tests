#!/usr/bin/env R

# Author: Sean Maden
#
# Format results output from IRFinder-S

srrid <- "SRR6026510"; run.handle <- "ipsc"

#----------
# load data
#----------
irfs.fname <- "IRFinder-IR-nondir.txt"
irfs <- data.table::fread(irfs.fname, header = T, data.table = F)

table(irfs$Warnings)
#
#     -              LowCover           LowSplicing          MinorIsoform NonUniformIntronCover 
# 66592                125971                  2623                 12835                 41958 

table(irfs$Warnings)/nrow(irfs)

#         -              LowCover           LowSplicing          MinorIsoform NonUniformIntronCover 
# 0.26639038            0.50392633            0.01049288            0.05134431            0.16784610

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