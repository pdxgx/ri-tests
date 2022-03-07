#!/usr/bin/env R

# Author: Sean Maden
# 
# Harmonize short and long read datasets

library(data.table)
library(GenomicRanges)

#----------
# load data
#----------
# long read data
tsv.fname <- "RI_txs_to_read_ids_final_09-21-2021_16.53.14.tsv"
tsv <- fread(tsv.fname, sep = "\t", header = T, data.table = F)

# load and format iread results
# load iread table
ir.fname <- "SRR6026510_allchr_iread.txt"
ir <- read.table(ir.fname, sep = " ", header = T)
# retain highest count among duplicated regions
dupv <- unique(ir[duplicated(ir[,1]),1])
irf <- ir[!ir[,1] %in% dupv,]
irdup <- ir[ir$intron_id %in% dupv,]
irdup <- irdup[order(irdup$intron_id, rev(irdup$fpkm)),]
irf <- rbind(irf, irdup[!duplicated(irdup[,1]),])
# save new table
irf.fname <- "SRR6026510_allchr_iread-filtdup"
save(irf, file = paste0(irf.fname, ".rda"))
write.table(irf, file = paste0(irf.fname, ".txt"))
write.csv(irf, file = paste0(irf.fname, ".csv"))

# load and format interest results
it.fname <- "SRR6026510_allchr_interest.tsv"
it <- fread(it.fname, header = F, data.table = F); it <- it[,c(2:ncol(it))]
colnames(it) <- c("chr", "start", "end", "strand",
                  "int_ex", "int_ex_num", "collapsed_transcripts_id",
                  "collapsed_transcripts", "collapsed_gene_id",
                  "IntRet_frequency", "IntRet_genewide_scaled")
nrow(it) # 617212
ftid <- paste0(it[,1], ":", it[,2], "-", it[,3])
length(unique(ftid)) # 617212
# save it
it.fname <- "SRR6026510_allchr_interest-format"
write.csv(it, file = paste0(it.fname, ".csv"))
write.table(it, file = paste0(it.fname, ".tsv"))
save(it, file = paste0(it.fname, ".rda"))

# load superintronic data
si.fname <- "superintronic_allchr_SRR6026510.rda"
si <- get(load(si.fname))
dim(si) # [1] 47282793       12

# kma/express
kma.fname <- "express_counts_SRR6026510.rda"
kma <- get(load(kma.fname))
kma.tpm <- kma$tpm
idv <- kma.tpm$target_id
kma.tpm.df <- data.frame(chr = gsub(":.*", "", idv),
                         start = gsub("-.*|.*:", "", idv),
                         end = gsub(".*-", "", idv),
                         kma.tpm = kma.tpm[,1], stringsAsFactors = F)
dim(kma.tpm.df) # [1] 438216      4
kma.tpm.df <- kma.tpm.df[!is.na(as.numeric(kma.tpm.df$start)),]
dim(kma.tpm.df) # [1] 208636      4

#-----------------
# helper functions
#-----------------
# calculate data summaries across intv.gr ranges
# using median values of col cname from data.df
grstat_all <- function(data.gr, cname = "intron_expression", 
                       new_cname = "intron_expression",
                       intv.gr = intv.gr){
  message("Finding overlaps between data.gr and intv.gr...")
  fol.gr <- GenomicRanges::findOverlaps(data.gr, intv.gr)
  sh <- subjectHits(fol.gr); unique.sh <- unique(sh)
  message("Summarizing variable ", cname, " medians on overlaps...")
  med.scorev <- unlist(lapply(unique.sh, function(shi){
    which.shi <- which(sh == shi)
    qhif <- unique(queryHits(fol.gr[which.shi]))
    median(data.gr@elementMetadata[,cname][qhif])}))
  data.intv.gr <- intv.gr[unique.sh]
  eval(parse(text = paste0("data.intv.gr@elementMetadata$",
                           new_cname," <- med.scorev")))
  return(data.intv.gr)
}





#-------------
# append iread
#-------------



# parse overlaps
fol.gr <- findOverlaps(ir.gr, intv.gr)
sh <- subjectHits(fol.gr); unique.sh <- unique(sh)
med.scorev <- unlist(lapply(unique.sh, function(shi){
  which.shi <- which(sh == shi)
  qhif <- unique(queryHits(fol.gr[which.shi]))
  median(ir.gr$fpkm[qhif])}))
ir.intv.gr <- intv.gr[unique.sh]
ir.intv.gr$ir_median_score <- med.scorev
ir.intv <- ir.intv.gr

#----------------
# append interest
#----------------




#--------------------
# get data as granges
#--------------------
# long read intron granges
intv <- unlist(lapply(unique(tsv$transcript), function(tid){
  message(tid); tsvf <- tsv[which(tsv$transcript == tid)[1],]
  intv <- unique(unlist(strsplit(tsvf$tx_introns, ";")))
  intv <- paste0(tsvf$chrom[1], ":", intv); return(intv)}))
df.intv <- data.frame(chr = gsub(":.*", "", intv),
                      start = gsub(".*:|-.*", "", intv),
                      end = gsub(".*-", "", intv), 
                      stringsAsFactors = F)
intv.gr <- makeGRangesFromDataFrame(df.intv)

# iread
# get the granges and find overlaps
ir.df <- irf; intronv <- ir.df$intron_id; lsplit <- strsplit(intronv, "-")
ir.df$chr <- paste0("chr", unlist(lapply(lsplit, 
                                         function(li){li[[1]]})))
ir.df$start <- unlist(lapply(lsplit, function(li){li[[2]]}))
ir.df$end <- unlist(lapply(lsplit, function(li){li[[3]]}))
ir.gr <- makeGRangesFromDataFrame(ir.df, keep.extra.columns = T)

# interest
# get the granges and find overlaps
it.intron <- it[it$int_ex == "intron",]
it.gr <- makeGRangesFromDataFrame(it.intron, keep.extra.columns = T)

# kma/express


# supreintronic
# get granges -- superintronic
num.rows <- 10000000; si.intv <- GRanges()
for(start.index in seq(1, nrow(si), num.rows)){
  end.index <- (start.index+num.rows-1)
  end.index <- ifelse(end.index > nrow(si), nrow(si), end.index)
  si.gr <- makeGRangesFromDataFrame(si[start.index:end.index,], 
                                    keep.extra.columns = T)
  fol.gr <- findOverlaps(si.gr, intv.gr)
  sh <- subjectHits(fol.gr); unique.sh <- unique(sh)
  med.scorev <- unlist(lapply(unique.sh, function(shi){
    which.shi <- which(sh == shi)
    qhif <- unique(queryHits(fol.gr[which.shi]))
    median(si.gr$score[qhif])}))
  si.intv.gr <- intv.gr[unique.sh]
  si.intv.gr$si_median_score <- med.scorev
  si.intv <- unlist(list(si.intv, si.intv.gr))
  message(start.index)}
si.intv.all <- do.call("c", si.intv) # bind ranges list
# save table
si.new.fname <- "SRR6026510_allchr_superintronic-byintron"
save(si.intv.all, file = paste0(si.new.fname, ".rda"))
write.table(si.intv.all, file = paste0(si.new.fname, ".txt"))
write.csv(si.intv.all, file = paste0(si.new.fname, ".csv"))

#------------------------------
# bind data to long read ranges
#------------------------------

# bind interest counts to long read ranges

# interest
it.intv <- grstat_all(it.intron, cname = "IntRet_frequency")

#-------------------
# append kma/express 
#-------------------


# make granages and bind
kma.tpm.gr <- makeGRangesFromDataFrame(kma.tpm.df, 
                                       keep.extra.columns = T)
length(kma.tpm.gr) # 208636
# parse overlaps
cname <- "kma.tpm"
fol.gr <- findOverlaps(kma.tpm.gr, intv.gr)
sh <- subjectHits(fol.gr); unique.sh <- unique(sh)
med.scorev <- unlist(lapply(unique.sh, function(shi){
  which.shi <- which(sh == shi)
  qhif <- unique(queryHits(fol.gr[which.shi]))
  median(kma.tpm.gr$kma.tpm[qhif])}))
kma.intv.gr <- intv.gr[unique.sh]
kma.intv.gr$kma_median_tpm <- med.scorev
kma.intv <- kma.intv.gr
save(kma.intv, file = "kma_intv_gr.rda")
# kma.intv <- grstat_all(kma.tpm.df, cname = "kma.tpm")

#-------------
# bind results
#-------------



