#!/usr/bin/env R

# Author: Sean Maden
#
# Pairs plots of each SR tool by sample.

library(GGally)

#----------
# load data
#----------
# get background table
plot.titlev <- c("HX1", "iPSC")
# tsv.fname.ipsc <- "called_RI_data_summary_iPSC.tsv"
# tsv.fname.hx1 <- "called_RI_data_summary_HX1.tsv"
tsv.fname.ipsc <- "called_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv"
tsv.fname.hx1 <- "called_RI_data_summary_HX1featureannotated_GCcontent.tsv"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

#---------------------------------
# pairs plots -- no transformation
#---------------------------------
# get plot objects
lgg <- lapply(seq(length(ltsv)), function(ii){
  tsvi <- ltsv[[ii]]; samplei <- names(ltsv)[ii]
  tsvi <- tsvi[!duplicated(tsvi$intron),]
  tsvi <- tsvi[,grepl(".*lwm$", colnames(tsvi))]
  colnames(tsvi) <- c("iREAD", "IntEREst", "superintronic", "KMA", "IRFinder-S",
                      "MAJIQ", "rMATS", "SUPPA2")
  colnames(tsvi) <- gsub("_.*", "", colnames(tsvi))
  ggpairs(tsvi, title= paste0(samplei),
          upper = list(continuous = wrap("points", alpha = 0.3), 
                       combo = wrap("dot_no_facet", alpha = 0.4)),
          lower = list(continuous = wrap("points", alpha = 0.3), 
                       combo = wrap("dot_no_facet", alpha = 0.4))) + 
    theme_bw()
}); names(lgg) <- names(ltsv)

# save figures
plot.fname.ipsc <- "ggpairs-srtool-cor_ipsc"
plot.fname.hx1 <- "ggpairs-srtool-cor_hx1"
width <- 15; height <- 15
# pdf
pdf(paste0(plot.fname.hx1, ".pdf"), width, height)
lgg[["HX1"]]; dev.off()
pdf(paste0(plot.fname.ipsc, ".pdf"), width, height)
lgg[["iPSC"]]; dev.off()
# png
# png(paste0(plot.fname.hx1, ".png"), width = width, height = height, 
#    units = "in", res = 500)
# lgg[["HX1"]]; dev.off()
# png(paste0(plot.fname.ipsc, ".png"), width = width, height = height, 
#    units = "in", res = 500)
# lgg[["iPSC"]]; dev.off()
