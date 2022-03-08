#!/usr/bin/env R

# Author: Sean Maden
# Format the kma/express results

library(GenomicRanges)
library(dplyr)

srrid <- "SRR6026510"; run.handle <- "ipsc"

#-------------
# save all kma
#-------------
kma.fname <- paste0("express_counts_",srrid,".rda")
kma <- get(load(kma.fname))

# get the formatted df
kma.tpm <- kma$tpm; idv <- kma.tpm$target_id
kma.tpm.df <- data.frame(chr = gsub(":.*", "", idv),
                         start = gsub("-.*|.*:", "", idv),
                         end = gsub(".*-", "", idv),
                         kma.tpm = kma.tpm[,1], stringsAsFactors = F)
dim(kma.tpm.df) 
# [1] 438216      4
kma.tpm.df <- kma.tpm.df[!is.na(as.numeric(kma.tpm.df$start)),]
dim(kma.tpm.df) 
# [1] 208636      4

# check for duplicates
kma.tpm.df$intronid <- paste0(kma.tpm.df$chr, ":", 
                              kma.tpm.df$start, "-", 
                              kma.tpm.df$end)
length(which(duplicated(kma.tpm.df$intronid))) 
# 0

#----------------------
# make granges and save
#----------------------
kma.gr <- makeGRangesFromDataFrame(kma.tpm.df, 
                                   keep.extra.columns = T)

# save granges
kma.gr.fname <- paste0("granges_kma_",srrid,"-",run.handle,".rda")
save(kma.gr, file = kma.gr.fname)

# save df
kma.gr.df.fname <- paste0("df-granges_kma_",
                          srrid,"-",run.handle,".csv")
write.csv(kma.tpm.df, file = kma.gr.df.fname)

#---------------------------------------------
# process kma intron retention ("ir") restults
#---------------------------------------------
# load kma/express *only retained introns*
kma.ri.fname <- paste0("ir-results-kma_",
                       srrid,"-",run.handle,".rda")
kma.ri <- get(load(kma.ri.fname))
kma.ri <- kma.ri$flat
# replace NaN retention with 0
which.retention <- is.na(kma.ri$retention)
kma.ri$retention[which.retention] <- 0
# check duplicated introns
length(which(duplicated(kma.ri$intron)))
# 0

# save the flat matrix
kma.ri.fname <- paste0("df-granges_kma-ir_",srrid,"-", run.handle, ".csv")
write.csv(kma.ri, file = kma.ri.fname)

# get granges
kma.ri$chr <- gsub(":.*", "", kma.ri$intron)
kma.ri$start <- gsub("(.*:|-.*)", "", kma.ri$intron)
kma.ri$end <- gsub(".*-", "", kma.ri$intron)
kma.ri.gr <- makeGRangesFromDataFrame(kma.ri, keep.extra.columns = T)
kma.ri.gr.fname <- paste0("granges_kma-ir_",srrid,"-",run.handle,".rda")
save(kma.ri.gr, file = kma.ri.gr.fname)
