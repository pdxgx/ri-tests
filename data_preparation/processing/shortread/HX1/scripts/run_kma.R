#!/usr/bin/env R

# Author: Sean Maden
#
# Script to process KMA outputs in R.
#
# Remote server instructions:
# 1. run R with sudo access
# > sudo /lib64/R/bin/$

library(kma)

# set run id
srrid <- "SRR2911306"

#----------
# load data
#----------
save.dpath <- file.path("home", "metamaden", "ri_results", srrid)
resources.dpath <- file.path("eternity", "data", "RI_benchmarking_results")
xprs.dpath <- file.path(resources.dpath, paste0("express_output_",srrid))
xprs.fnv <- list.files(xprs.dpath)
xprs.results.fpath <- file.path(xprs.dpath, "results.xprs")
xprs <- read_express(xprs.results.fpath, srrid, "ipsc")
# load intron_to_trans
it.fpath <- file.path(resources.dpath, "gencode_v35_annotation_files", 
                      "sm", "intron_to_transcripts.txt")
intron_to_trans <- data.table::fread(it.fpath, data.table = FALSE)

#-------------------
# save results table
#-------------------
# save express robject
xp.fname <- paste0("express_counts_",srrid,".rda")
xp.fpath <- file.path(save.dpath, save.fname)
save(xprs, file = xp.fpath)

#--------------------------
# quantify intron retention
#--------------------------
# make new intron retention object
ir.fname <- paste0("ir-results-kma_",srrid,".rda")
ir.fpath <- file.path(save.dpath, ir.fname)
ir <- newIntronRetention(xprs$tpm, intron_to_trans, 
                         xprs$condition, xprs$uniq_counts)
# save intron retention results
save(ir, file = ir.fpath)