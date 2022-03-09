#!/usr/bin/env R

# Author: Sean Maden
# 
# Harmonize short and long read datasets. This script calculates the length-weighted
# medians (LWMs) across introns with quantified expression. Note that LWMs are 
# to introns expressed in the long read samples. The script further identifies the 
# subset of expressed introns meeting tool-specific filter criteria, a.k.a. "
# called RIs".
# 

library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# set run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

#----------
# load data
#----------
# long read granges
lr.gr.fname <- paste0("granges_longread_",srrid,"-",run.handle,".rda")
intv.gr <- get(load(lr.gr.fname))

#-----------------
# helper functions
#-----------------
# get length weighted medium
get_lwm <- function(val, weights){
  if(!length(val) == length(weights)){
    stop("Error: length of val must equal length of weights.")}
  return(median(rep(val, times = weights)))
}

# bind new length weighted mediums
grbind_summary <- function(data.gr, intv.gr = intv.gr, cname = "intron_expression", 
                           new.cname = "intron_expression"){
  if(!cname %in% colnames(data.gr@elementMetadata)){
    stop("Provided variable name ",cname," not in data.gr@elementMetadata")}
  message("Finding overlaps between data.gr and intv.gr...")
  fol.gr <- findOverlaps(data.gr, intv.gr)
  sh <- subjectHits(fol.gr); qh <- queryHits(fol.gr)
  valv <- data.gr@elementMetadata[,cname][qh]
  dataf.gr <- data.gr[qh]
  dfol <- data.frame(width = abs(start(dataf.gr)-end(dataf.gr)), 
                     val = valv, stringsAsFactors = F)
  dfol <- dfol[dfol$val > 0,]
  dfol <- data.frame(sh = sh, qh = qh, val = valv, stringsAsFactors = F)
  dfol$width <- width(data.gr)[qh] # append region widths
  message("Getting the region value summaries...")
  dcol <- dfol %>% group_by(sh) %>% # collapse on subject hits
    summarise(lwm = get_lwm(val, weights = width))
  new.dat.wmedian <- rep(0, length(intv.gr))
  new.dat.wmedian[dcol$sh] <- dcol$lwm
  txt.wmedian <- paste0("intv.gr@elementMetadata$",new.cname,"_lwm <- new.dat.wmedian")
  eval(parse(text = txt.wmedian))
  return(intv.gr)
}

#------------------------------------------------------
# map SR to LR ranges, all introns and filtered introns
#------------------------------------------------------
# iREAD
# Notes: Filter on fragments >= 20, junction reads >= 0.9, 
# FPKM >= 3, and entropy score >= 0.9 (per tool defaults/paper)
ir.gr <- get(load(paste0("granges_iread_",srrid,"-",run.handle,".rda")))
# map all ranges
intv.gr <- grbind_summary(data.gr = ir.gr, intv.gr = intv.gr, cname = "fpkm", 
                          new.cname = "iread_fpkm_allintron")
# map filtered ranges
which.ir.gr <- which(ir.gr$fragments >= 10 & ir.gr$junction_reads >= 1 &
                       ir.gr$fpkm >= 1 & ir.gr$entropy_score >= 0.9)
intv.gr <- grbind_summary(data.gr = ir.gr[which.ir.gr], intv.gr = intv.gr, 
                          cname = "fpkm", new.cname = "iread_fpkm_filtintron")

# IntEREst
# Notes: Filter on FPKM >= 45 (per target order of magnitude)
it.gr <- get(load(paste0("granges_interest_",srrid,"-",run.handle,".rda")))
# map all ranges
intv.gr <- grbind_summary(data.gr = it.gr, intv.gr = intv.gr, cname = "IntRet_frequency", 
                          new.cname = "interest_fpkm_allintron")
# map filtered ranges
intv.gr <- grbind_summary(data.gr = it.gr[it.gr$IntRet_frequency >= 45], intv.gr = intv.gr, 
                          cname = "IntRet_frequency", new.cname = "interest_fpkm_filtintron")

# superintronic
# Notes: Filter on LWM >= 5 by LR region
si.gr <- get(load(paste0("granges_superintronic_",srrid,"-",run.handle,".rda")))
si.gr <- si.gr[si.gr$feature_type=="intron"]
# map all ranges
intv.gr <- grbind_summary(data.gr = si.gr, intv.gr = intv.gr, cname = "score", 
                          new.cname = "superintronic_score_allintron")
# map filtered ranges
si.grf <- subsetByOverlaps(si.gr, intv.gr[intv.gr$superintronic_score_allintron_lwm >= 5])
intv.gr <- grbind_summary(data.gr = si.grf, intv.gr = intv.gr, cname = "score", 
                          new.cname = "superintronic_score_filtintron")
rm(si.gr) # remove large in-mem object

# KMA
# Notes: Filter on intron read counts >= 10
kma.gr <- get(load(paste0("granges_kma-express_",srrid,".rda")))
kma.ir.gr <- get(load(paste0("granges_kma-ir_", srrid, "-", run.handle, ".rda")))
# map all ranges
intv.gr <- grbind_summary(data.gr = kma.gr, intv.gr = intv.gr, 
                          cname = "kma.tpm", new.cname = "kma_tpm_allintron")
# map filtered ranges
kma.grf <- subsetByOverlaps(kma.gr, kma.ir.gr[kma.ir.gr$unique_counts >= 10])
intv.gr <- grbind_summary(data.gr = kma.grf, intv.gr = intv.gr, 
                          cname = "kma.tpm", new.cname = "kma_tpm_filtintron")

# IRFinder-S
# Notes: Filter on expression flag and IRratio > 0.05
irfs.gr <- get(load(paste0("granges_irfinders_", srrid,"-",run.handle,".rda")))
# map all ranges
intv.gr <- grbind_summary(data.gr = irfs.gr, intv.gr = intv.gr, cname = "IRratio", 
                          new.cname = "irfinders_irratio_allintron")
# map filtered ranges
irfs.grf <- irfs.gr[irfs.gr$Warnings=="-" & irfs.gr$IRratio > 0.05]
intv.gr <- grbind_summary(data.gr = irfs.grf, intv.gr = intv.gr, cname = "IRratio", 
                          new.cname = "irfinders_irratio_filtintron")

#-------------
# save results
#-------------
# write new table
csv.fname <- paste0("granges-lrmap_sr-5-methods_", srrid,"-",run.handle,".csv")
write.csv(as.data.frame(intv.gr), file = csv.fname, row.names = F)

# save new bound data
gr.fname <- paste0("granges-lrmap_sr-5-methods_", srrid,"-",run.handle,".rda")
save(intv.gr, file = gr.fname)

#--------------------------------------
# barplot summaries, all vs. called RIs
#--------------------------------------
md <- intv.gr@elementMetadata; mdf <- md[,grepl("lwm$", colnames(md))]
toolv <- unique(gsub("_.*", "", colnames(mdf))); cnv <- colnames(mdf)
dfp <- do.call(rbind, lapply(toolv, function(tooli){
  do.call(rbind, lapply(c("allintron", "filtintron"), function(lvli){
    which.cni <- cnv[grepl(tooli, cnv) & grepl(lvli, cnv)]
    numi <- length(which(md[,which.cni] > 0))
    data.frame(tool = tooli, type = lvli, num.introns = numi)
  }))
}))

# format vars
dfp$level <- paste0(dfp$tool, ";", dfp$type)
dfp$Tool <- ifelse(dfp$tool == "interest", "IntEREst",
                   ifelse(dfp$tool == "kma", "KMA",
                          ifelse(dfp$tool == "iread", "iREAD",
                                 ifelse(dfp$tool == "irfinders", "IRFinder-S",
                                        dfp$tool))))

# color palette
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7')

# get plot object
bp <- ggplot(dfp, aes(x = level, y = num.introns, fill = Tool)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  geom_text(aes(x = level, y = num.introns + 1000, label = num.introns)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1)) +
  scale_fill_manual(values = pal) + ggtitle("iPSC")

# save new figure
pdf.fname <- paste0("bp_intron-counts_bytool_", srrid,"-",run.handle,".pdf")
pdf(pdf.fname, 6.5, 4.5); print(bp); dev.off()
