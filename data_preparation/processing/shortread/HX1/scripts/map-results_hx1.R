#!/usr/bin/env R

# Author: Sean Maden
# 
# Harmonize short and long read datasets

library(data.table)
library(GenomicRanges)
library(dplyr)

srrid <- "SRR2911306"

#-----------------
# helper functions
#-----------------
grbind_summary <- function(data.gr, intv.gr = intv.gr, 
                          cname = "intron_expression", 
                          new.cname = "intron_expression"){
  if(!cname %in% colnames(data.gr@elementMetadata)){
    stop("Provided variable name ",cname," not in data.gr@elementMetadata")}
  message("Finding overlaps between data.gr and intv.gr...")
  fol.gr <- findOverlaps(data.gr, intv.gr)
  sh <- subjectHits(fol.gr); qh <- queryHits(fol.gr)
  valv <- data.gr@elementMetadata[,cname][qh]
  dfol <- data.frame(sh = sh, qh = qh, val = valv, stringsAsFactors = F)
  message("Getting the region value summaries...")
  dcol <- dfol %>% group_by(sh) %>% # collapse on subject hits
    summarise(median = median(val, na.rm = TRUE),
              mean = mean(val, na.rm = TRUE))
  new.dat.median <- new.dat.mean <- rep(0, length(intv.gr))
  new.dat.median[dcol$sh] <- dcol$median; new.dat.mean[dcol$sh] <- dcol$mean
  txt.median <- paste0("intv.gr@elementMetadata$",new.cname,"_median <- new.dat.median")
  txt.mean <- paste0("intv.gr@elementMetadata$",new.cname,"_mean <- new.dat.mean")
  eval(parse(text = txt.median)); eval(parse(text = txt.mean))
  return(intv.gr)
}

#-----------------------------------------
# load granges and map to long read ranges
#-----------------------------------------
# long read granges
lr.gr.fname <- "granges_longread_hx1.rda"
intv.gr <- get(load(lr.gr.fname))
intv.gr <- reduce(intv.gr)

# load formatted iread results
ir.gr.fname <- paste0("granges_iread_",srrid,".rda")
ir.gr <- get(load(ir.gr.fname))
intv.gr <- grbind_summary(data.gr = ir.gr, intv.gr = intv.gr, 
                     cname = "fpkm", new.cname = "iread_fpkm")

# load formatted interest results
it.gr.fname <- paste0("granges_interest_",srrid,".rda")
it.gr <- get(load(it.gr.fname))
intv.gr <- grbind_summary(data.gr = it.gr, intv.gr = intv.gr, 
                     cname = "IntRet_frequency", 
                     new.cname = "interest_intret_freq")

# load superintronic granges
si.gr.fname <- paste0("granges_superintronic_",srrid,".rda")
si.gr <- get(load(si.gr.fname))
intv.gr <- grbind_summary(data.gr = si.gr, intv.gr = intv.gr, 
                     cname = "si_median_score", 
                     new.cname = "si_median_score")

# kma/express
kma.gr.fname <- paste0("granges_kma-express_",srrid,".rda")
kma.gr <- get(load(kma.gr.fname))
intv.gr <- grbind_summary(data.gr = kma.gr, intv.gr = intv.gr, 
                      cname = "kma.tpm", new.cname = "kma_tpm")
#-----
# save
#-----
# save new bound data
gr.fname <- paste0("granges-longread_score-mean-median_ir-4-methods_",srrid,"-hx1.rda")
save(intv.gr, file = gr.fname)

# save csv
gr.csv <- as.data.frame(intv.gr, stringsAsFactors = F)
gr.csv.fname <- paste0("csv_longread-ranges_ir-4-methods_",srrid,"-hx1.csv")
write.csv(gr.csv, file = gr.csv.fname)

#-----
# cor
#----
mcor.median <- gr.csv[,grepl("median$", colnames(gr.csv))]
mcor.mean <- gr.csv[,grepl("mean$", colnames(gr.csv))]
cor(mcor.median, method = "spearman")
# iread_fpkm_median interest_intret_freq_median si_median_score_median
# iread_fpkm_median                  1.00000000                          NA             0.01537794
# interest_intret_freq_median                NA                           1                     NA
# si_median_score_median             0.01537794                          NA             1.00000000
# si_median_score_mean               0.01569570                          NA             0.99620829
# kma_tpm_median                     0.14234819                          NA             0.18791300
# si_median_score_mean kma_tpm_median
# iread_fpkm_median                      0.0156957      0.1423482
# interest_intret_freq_median                   NA             NA
# si_median_score_median                 0.9962083      0.1879130
# si_median_score_mean                   1.0000000      0.1874976
# kma_tpm_median                         0.1874976      1.0000000
cor(mcor.mean, method = "spearman")
# iread_fpkm_mean interest_intret_freq_mean si_median_score_mean
# iread_fpkm_mean                1.00000000                 0.1187795           0.01585838
# interest_intret_freq_mean      0.11877951                 1.0000000           0.39914006
# si_median_score_mean           0.01585838                 0.3991401           1.00000000
# kma_tpm_mean                   0.14245636                 0.6005624           0.18872193
# kma_tpm_mean
# iread_fpkm_mean              0.1424564
# interest_intret_freq_mean    0.6005624
# si_median_score_mean         0.1887219
# kma_tpm_mean                 1.0000000
