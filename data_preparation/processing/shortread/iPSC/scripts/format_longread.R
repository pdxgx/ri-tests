#!/usr/bin/env R

# Author: Sean Maden
# 
# Format the long read intron ranges table.
#

library(data.table)

# set run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

#----------
# load data
#----------
# get the long read ranges
tsv.fname <- "RI_txs_to_read_ids_final.tsv"
tsv <- fread(tsv.fname, sep = "\t", header = T, data.table = F)

#---------------------------
# get data for intron ranges
#---------------------------
tsv$read_left <- as.numeric(tsv$read_left)
tsv$read_right <- as.numeric(tsv$read_right)
tidv <- unique(tsv$transcript)
dfexpr <- do.call(rbind, lapply(tidv, function(tid){
  tf <- tsv[tsv$transcript == tid,]; chr <- unique(tf$chrom)
  intv <- unlist(strsplit(unique(tf$tx_introns), ";"))
  intv.str <- paste0(chr, ":", intv);fract.denom <- nrow(tf)
  intv.fractexpr.numreads <- do.call(rbind, lapply(intv, function(inti){
    int.start <- as.numeric(gsub("-.*", "", inti))
    int.end <- as.numeric(gsub(".*-", "", inti))
    which.coord <- tf$read_left <= int.start & tf$read_right >= int.end
    tff <- tf[which.coord,,drop = F]; expri <- "NA"; num.reads <- 0
    if(nrow(tff) > 0){
      num.reads <- nrow(tff); inti.str <- paste0("(^|;)", inti, "($|;)")
      # get fraction of reads having the intron 
      # e.g. for which inti is not listed under read_introns
      which.inti.expr <- !grepl(inti.str, tff$read_introns)
      fract.numi <- nrow(tff[which.inti.expr,])
      expri <- fract.numi/fract.denom}; return(c(expri, num.reads))}))
  data.frame(tid = rep(tid, length(intv)), intid = intv.str,
             num.reads = intv.fractexpr.numreads[,2],
             fract.expr = intv.fractexpr.numreads[,1], 
             stringsAsFactors = F)}))
# format vars
for(c in 3:4){dfexpr[,c] <- as.numeric(dfexpr[,c])}

# summarize by intron id
dfint <- dfexpr %>% group_by(intid) %>% # collapse on subject hits
  summarise(fract.expr.median = median(fract.expr, na.rm = TRUE),
            fract.expr.mean = mean(fract.expr, na.rm = TRUE),
            min.expr = min(fract.expr, na.rm = TRUE),
            max.expr = max(fract.expr, na.rm = TRUE),
            numreads.median = median(num.reads, na.rm = TRUE),
            numreads.mean = mean(num.reads, na.rm = TRUE))
dfint <- as.data.frame(dfint, stringsAsFactors = F)

#---------------------
# get granges and save
#---------------------
# get the intron granges
intv <- as.character(dfint[,1])
dfgr <- data.frame(chr = gsub(":.*", "", intv),
                   start = gsub(".*:|-.*", "", intv),
                   end = gsub(".*-", "", intv), 
                   stringsAsFactors = F)
dfgr <- cbind(dfgr, dfint[,c(2:7)])
intv.gr <- makeGRangesFromDataFrame(dfgr, 
                                    keep.extra.columns = T)

# save granges
intv.gr.fname <- paste0("granges_longread_", 
                        srrid,"-", run.handle, ".rda")
save(intv.gr, file = intv.gr.fname)