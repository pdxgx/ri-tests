#!/usr/bin/env R

# Author: Sean Maden
# Format the IntEREst results.

library(dplyr)
library(GenomicRanges)

srrid <- "SRR2911306"; run.handle <- "hx1"

# load the results table
it.fname <- paste0("interest_allchr_",srrid,".rda")
it <- get(load(it.fname))
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
length(dup)
# [1] 7015680

# get new df for duplicated ids
which.notdup <- it.df$intronid %in% dup.id & !duplicated(it.df$intronid)
df.notdup <- it.df[which.notdup,]
dim(df.notdup) 
# [1] 292320     12
# retain max observed expr by dup id
dfply <- it.df[it.df$intronid %in% dup.id,]; 
dim(dfply) 
# [1] 7308000      12
eval.str <- paste0("dfply <- dfply %>% group_by(intronid) ",
                   "%>% summarise(",cname," = max(",cname,", na.rm = T))")
eval(parse(text = eval.str))
dim(dfply) 
# [1] 292320      2
dfply <- dfply[order(match(dfply$intronid, df.notdup$intronid)),]
identical(dfply$intronid, df.notdup$intronid) 
# TRUE
df.notdup[,cname] <- dfply[,2]
# bind results
it.dff <- rbind(it.df[!it.df$intronid %in% dup.id,], df.notdup)
dim(it.df) 
# [1] 7308000      12
dim(it.dff) 
# [1] 292320     12

#---------------------
# get granges and save
#---------------------
it <- it.dff
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
