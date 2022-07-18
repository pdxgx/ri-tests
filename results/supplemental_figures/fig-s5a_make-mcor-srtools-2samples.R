#!/usr/bin/env R

# Author: Sean Maden
# 
# Analyze the short read data, mapped to long read intron ranges.

library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(GenomicRanges)

#----------
# load data
#----------
srridv <- c("SRR2911306", "SRR6026510")
run.handlev <- c("hx1", "ipsc")
lgrpile <- lirv <- list()
lsamp <- lapply(seq(length(srridv)), function(ii){
  srrid <- srridv[ii]; run.handle <- run.handlev[ii]
  dname <- paste0(srrid, "_", run.handle)
  # irv.fname <- paste0("granges-lrmap_sr-5-methods_", srrid,"-",run.handle,".rda")
  irv.fname <- paste0("granges-lrmap_sr-8-methods_", srrid,"-",run.handle,".rda")
  irv <- get(load(file.path(dname, irv.fname)))
  # get median gene coverage
  mpile.fname <- paste0("mpile-sstat_gencode-v35_",srrid,"-",run.handle,".csv")
  mpile <- read.csv(file.path(dname, mpile.fname))
  grpile <- makeGRangesFromDataFrame(mpile, keep.extra.columns = T)
  grpile <- grpile[!is.na(grpile$median)]
  grpile <- grpile[grpile$median >= 2]
  irvf <- subsetByOverlaps(irv, grpile) # filt irv on coverage
  list(grpile = grpile, irv = irv, irvf = irvf)
}); names(lsamp) <- run.handlev

# get intron data
plot.titlev <- c("HX1", "iPSC")
# tsv.fname.hx1 <- "target_genes_LR_annotated_granges-lrmap_sr-5-methods_SRR2911306-hx1.csv"
# tsv.fname.ipsc <- "target_genes_LR_annotated_granges-lrmap_sr-5-methods_SRR6026510-ipsc.csv"
tsv.fname.hx1 <- "target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR2911306-hx1.csv"
tsv.fname.ipsc <- "target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR6026510-ipsc.csv"
ltsv <- list()
ltsv[["iPSC"]] <- read.csv(tsv.fname.ipsc, header = T)
ltsv[["HX1"]] <- read.csv(tsv.fname.hx1, header = T)

#--------------
# get plot data
#--------------
# get the correlation matrices
lmcor <- lapply(seq(length(ltsv)), function(ii){
  run.handle <- names(ltsv)[ii]; tsvi <- ltsv[[ii]]
  # get metadata df
  which.cnv <- grepl("lwm", colnames(tsvi)) & grepl("filtintron", colnames(tsvi))
  dati.cor <- tsvi[!duplicated(tsvi$intron),]
  dati.cor <- dati.cor[,which.cnv]
  dati.cor[is.na(dati.cor)] <- 0
  return(cor(dati.cor, method = "spearman"))
})
names(lmcor) <- tolower(names(ltsv))

# get dfp list
ldfp <- lapply(seq(length(lmcor)), function(ii){
  dfp <- lmcor[[ii]]
  rownames(dfp) <- colnames(dfp) <- c("iREAD", "IntEREst", "superintronic", 
                                      "KMA", "IRFinder-S", "MAJIQ", "rMATS", 
                                      "SUPPA2")
  if(names(lmcor)[ii] == "hx1"){
    which.na <- which(upper.tri(dfp)|dfp==1); dfp[which.na] <- "NA"
  } else{which.na <- which(lower.tri(dfp)|dfp==1); dfp[which.na] <- "NA"}
  dsh <- reshape2::melt(dfp); dsh$sample <- names(lmcor)[ii]
  dsh$Rho <- as.numeric(dsh$value)
  dsh$label <- as.character(round(dsh$Rho, 2)); dsh
}); names(ldfp) <- names(lmcor)

# get new plot object
dfp <- do.call(rbind, lapply(ldfp, function(dfi){dfi[!dfi$value == "NA",]}))
dfp$label <- as.character(round(as.numeric(dfp$value), 2))
dfp$Sample <- dfp$sample

# add manual color for hx1
# get the fill pattern legend
plot.legend1 <- ggplot(dfp, aes(x = Var1, y = Var2, fill = Rho)) + 
  theme_bw() + geom_tile() + 
  scale_fill_gradient2(low = "#520B71", mid = "white", 
                       high = "#E68400", midpoint = 0, 
                       name = expression(rho))
plot.legend1 <- get_legend(plot.legend1)
# get the ouline pattern legend
plot.legend2 <- ggplot(dfp, aes(x = Var1, y = Var2, color = Sample)) + geom_tile(size = 2, fill = "white") + 
  scale_color_manual(values = c("iPSC" = 'navyblue', "HX1" = "firebrick4"))
plot.legend2 <- get_legend(plot.legend2)

# get plot object
ggtile <- ggplot(dfp, aes(x = Var1, y = Var2, fill = Rho, color = sample)) +
  geom_tile(size = 1) + geom_text(aes(label = label), color = "black") +
  scale_color_manual(values = c("ipsc" = 'navyblue', "hx1" = "firebrick4")) +
  scale_fill_gradient2(low = "#520B71", mid = "white", high = "#E68400", midpoint = 0) +
  theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), 
                     axis.title.y = element_blank())

# save new plots# new pdf
plot.fname <- "mcor_srtools_2-samples"
pdf.fname <- paste0(plot.fname, ".pdf")
lm <- matrix(c(rep(1, 7),2,3), nrow = 1)
width <- 8; height <- 1.8
pdf(pdf.fname, width, height)
grid.arrange(ggtile, as_ggplot(plot.legend1), 
             as_ggplot(plot.legend2), layout_matrix = lm)
dev.off()
# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = width, height = height, units = "in", res = 500)
grid.arrange(ggtile, as_ggplot(plot.legend1), 
             as_ggplot(plot.legend2), layout_matrix = lm)
dev.off()

