#!/usr/bin/env R

# Author: Sean Maden
#
# Precision and recall by intron length

library(ggplot2)
library(ggpubr)
library(gridExtra)

#----------
# load data
#----------
plot.titlev <- c("HX1", "iPSC")
tsv.fname.hx1 <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR2911306-hx1.csv"
tsv.fname.ipsc <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR6026510-ipsc.csv"

ltsv <- list()
ltsv[["iPSC"]] <- read.csv(tsv.fname.ipsc)
ltsv[["HX1"]] <- read.csv(tsv.fname.hx1)

# get the color palette:
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')

#-----------------
# helper functions
#-----------------
# get precision and recall from a dfp of truth metrics
dfpr_bylenfilt <- function(tsv, len.filt, min.lr = 0.1, 
                           intron.type.cname = "filtintron", 
                           lr.metric.cname = "max_intron_persistence",
                           tool.strv = c("interest", "superintronic", "iread", 
                                         "kma", "irfinders", "rmats", "majiq",
                                         "suppa2")){
  do.call(rbind, lapply(len.filt, function(min.len){
    do.call(rbind, lapply(tool.strv, function(tooli){
      # get main df
      tsvi <- tsv[!duplicated(tsv$intron),]
      cnv <- colnames(tsvi)
      cnv <- c(cnv[grepl(tooli, cnv)], lr.metric.cname, "width")
      tsvi <- tsvi[,cnv]
      tsvi <- tsvi[tsvi$width >= min.len,] # filter on length
      which.tooli.value <- grepl(".*weighted_median$|.*lwm$", cnv)
      which.tooli.value <- which.tooli.value & grepl(intron.type.cname, cnv)
      tsvi[is.na(tsvi[,which.tooli.value]), which.tooli.value] <- 0
      # get prec/rec
      tp <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] > 0))
      fp <- length(which(tsvi[,lr.metric.cname] < min.lr & tsvi[,which.tooli.value] > 0))
      fn <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] == 0))
      precisioni <- ifelse((tp > 0 | fp > 0), tp/(tp+fp), 0)
      recalli <- ifelse((tp > 0 | fn > 0), tp/(tp+fn), 0)
      fm <- ifelse((precisioni > 0|recalli > 0), (2*precisioni*recalli)/(precisioni+recalli), 0)
      data.frame(tool = tooli, min.len = min.len,
                 precision = precisioni, recall = recalli, 
                 fmeasure = fm, num.lr = nrow(tsvi))
    }))
  }))
}

# get precision and recall from a dfp of truth metrics, with binning
dfpr_bylenfilt_binned <- function(tsv, len.filt, min.lr = 0.1, bin.size = 500,
                                  intron.type.cname = "filtintron", 
                                  lr.metric.cname = "max_intron_persistence",
                                  tool.strv = c("interest", "superintronic", 
                                                "iread", "kma", "irfinders",
                                                "rmats", "majiq", "suppa2")){
  do.call(rbind, lapply(seq(length(len.filt)), function(ii){
    do.call(rbind, lapply(tool.strv, function(tooli){
      # get main df
      tsvi <- tsv[!duplicated(tsv$intron),]
      cnv <- colnames(tsvi)
      cnv <- c(cnv[grepl(tooli, cnv)], lr.metric.cname, "width")
      # get len filtered data
      min.len <- len.filt[ii]; max.len <- min.len + bin.size
      which.sizefilt <- which(tsvi$width >= min.len & tsvi$width < max.len)
      tsvi <- tsvi[which.sizefilt,] # filter on length
      num.intron <- length(unique(tsvi$intron)) # get unique introns
      tsvi <- tsvi[,cnv] # get filtered columns
      # replace missing values
      which.tooli.value <- grepl(".*weighted_median$|.*lwm$", cnv)
      which.tooli.value <- which.tooli.value & grepl(intron.type.cname, cnv)
      tsvi[is.na(tsvi[,which.tooli.value]), which.tooli.value] <- 0
      # get prec/rec
      tp <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] > 0))
      fp <- length(which(tsvi[,lr.metric.cname] < min.lr & tsvi[,which.tooli.value] > 0))
      fn <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] == 0))
      precisioni <- ifelse((tp > 0 | fp > 0), tp/(tp+fp), 0)
      recalli <- ifelse((tp > 0 | fn > 0), tp/(tp+fn), 0)
      fm <- ifelse((precisioni > 0|recalli > 0), (2*precisioni*recalli)/(precisioni+recalli), 0)
      data.frame(tool = tooli, min.len = min.len, max.len = max.len, 
                 num.intron = num.intron, precision = precisioni, recall = recalli, 
                 fmeasure = fm, num.lr = nrow(tsvi))
    }))
  }))
}

#--------------
# get plot data
#--------------
# define the length filters
len.filt <- as.numeric(c(0, 
                         quantile(ltsv[[1]]$width, 
                                  seq(0, 0.9, 0.1), na.rm = T)))[1:10]
# get list of dfpr objects
ldfpr <- lapply(ltsv, function(tsvi){
  dfpr <- dfpr_bylenfilt(tsvi, len.filt = len.filt)
  # format vars
  dfpr[dfpr$tool=="interest",]$tool <- "IntEREst"
  #dfpr[dfpr$tool=="si",]$tool <- "superintronic"
  dfpr[dfpr$tool=="iread",]$tool <- "iREAD"
  dfpr[dfpr$tool=="kma",]$tool <- "KMA"
  dfpr[dfpr$tool=="irfinders",]$tool <- "IRFinder-S"
  dfpr$`RI detection\ntool` <- dfpr$tool; dfpr
})

# get the color palette:
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')
# set transparency
alpha.val <- 1
# set point and line size
pt.size <- 1.2; line.size <- 1
lgg <- lapply(ldfpr, function(dfpr){
  # format vars
  dfpr$Tool <- dfpr$tool
  # make plot objects
  gg.ptline.prec <- ggplot(dfpr, aes(x = min.len, y = precision, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab("Minimum intron length (bp)") +
    ylab("Precision") + scale_color_manual(values = pal)
  gg.ptline.rec <- ggplot(dfpr, aes(x = min.len, y = recall, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab("Minimum intron length (bp)") +
    ylab("Recall") + scale_color_manual(values = pal)
  gg.ptline.f1score <- ggplot(dfpr, aes(x = min.len, y = fmeasure, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab("Minimum intron length (bp)") +
    ylab("F1-score") + scale_color_manual(values = pal)
  plot.legend <- ggplot(dfpr, aes(x = min.len, y = recall, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + scale_color_manual(values = pal)
  plot.legend <- get_legend(plot.legend)
  # prep multiplot
  gg.ptline.prec <- gg.ptline.prec + theme(axis.title.x = element_blank(), legend.position = "none")
  gg.ptline.rec <- gg.ptline.rec + theme(axis.title.x = element_blank(), legend.position = "none")
  gg.ptline.f1score <- gg.ptline.f1score + theme(axis.title.x = element_blank(), legend.position = "none")
  return(list(prec = gg.ptline.prec, rec = gg.ptline.rec, 
              f1score = gg.ptline.f1score, legend = plot.legend))
})
names(lgg) <- names(ldfpr)

# plot dims
prec.ymax <- 1
rec.ymax <- 1
f1.ymax <- 1
xmax <- max(len.filt)

# get plot objects
prec.hx1 <- lgg[[1]]$prec + ylim(0, prec.ymax) + xlim(0, xmax) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
rec.hx1 <- lgg[[1]]$rec + ylim(0, rec.ymax) + xlim(0, xmax) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
f1.hx1 <- lgg[[1]]$f1score + ylim(0, f1.ymax) + xlim(0, xmax) +
  theme(axis.title.x = element_blank())
prec.ipsc <- lgg[[2]]$prec + ylim(0, prec.ymax) + xlim(0, xmax) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank())
rec.ipsc <- lgg[[2]]$rec + ylim(0, rec.ymax) + xlim(0, xmax) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank())
f1.ipsc <- lgg[[2]]$f1score + ylim(0, f1.ymax) + xlim(0, xmax) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank())

#-------------------------------------------------
# ggpt/line figure -- prec/rec/f1 by intron length
#-------------------------------------------------
# save new figure
# get figure vars
num.plot <- 4; num.legend <- 2
lm <- matrix(c(rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(2, num.plot + 2), rep(5, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend)),
             byrow = T, nrow = 3)
title.str <- paste0(paste0(rep(" ", 32), collapse = ""), 
                    "HX1", paste0(rep(" ", 52), collapse = ""), 
                    "iPSC", paste0(rep(" ", 90), collapse = ""),
                    collapse = "")
xlab.str <- paste0("Minimum intron length (bp)", paste0(rep(" ", 10), collapse = ""), collapse = "")
# save new plots
plot.fname <- "ggptline_prec-rec-by-ilength_combined-2samples"
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 7, 3.5)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm,
             top = title.str, bottom = xlab.str)
dev.off()
# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 7, height = 3.5, units = "in", res = 500)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm,
             top = title.str, bottom = xlab.str)
dev.off()

#------------------------------------------------
# get prec/rec/f1 for binned intron length ranges
#------------------------------------------------
# get data for fig 4b
# define the length filters
seq.max <- 4000; bin.size <- 250

len.filt <- seq(0, seq.max, bin.size)
# get list of dfpr objects
ldfpr <- lapply(ltsv, function(tsvi){
  dfpr <- dfpr_bylenfilt_binned(tsvi, len.filt = len.filt, bin.size = bin.size)
  # format vars
  dfpr[dfpr$tool=="interest",]$tool <- "IntEREst"
  dfpr[dfpr$tool=="iread",]$tool <- "iREAD"
  dfpr[dfpr$tool=="kma",]$tool <- "KMA"
  dfpr[dfpr$tool=="irfinders",]$tool <- "IRFinder-S"
  dfpr[dfpr$tool=="rmats",]$tool <- "rMATS"
  dfpr[dfpr$tool=="majiq",]$tool <- "MAJIQ"
  dfpr[dfpr$tool=="suppa2",]$tool <- "SUPPA2"
  dfpr$Tool <- dfpr$tool # format vars
  dfpr$`RI detection\ntool` <- dfpr$tool; dfpr
})
names(ldfpr) <- names(ltsv)

# set transparency
# alpha.val <- 1
# set point and line size
# pt.size <- 1.2; line.size <- 1

# plot params
# larger fonts
xaxis.text.size <- yaxis.text.size <- legend.text.size <- 14
title.text.size <- xaxis.title.size <- yaxis.title.size <- 18
# plot dims
prec.ymax <- 0.7
rec.ymax <- 1
f1.ymax <- 0.3
xmax <- max(len.filt) + 100
min.introns <- 80 # minimum introns filter term

lgg <- lapply(ldfpr, function(dfpr){
  dfprf <- dfpr[dfpr$num.intron >= min.introns,] # filter bins on min introns
  # make plot objects
  gg.prec <- ggplot(dfprf, aes(x = min.len, y = precision, color = Tool)) +
    geom_smooth(se = F) + theme_bw() + scale_color_manual(values = pal) +
    xlab("Intron length (# bp)") + ylab("Prec.")
  gg.rec <- ggplot(dfprf, aes(x = min.len, y = recall, color = Tool)) +
    geom_smooth(se = F) + theme_bw() + scale_color_manual(values = pal) +
    xlab("Intron length (# bp)") + ylab("Rec.")
  gg.f1 <- ggplot(dfprf, aes(x = min.len, y = fmeasure, color = Tool)) + 
    geom_smooth(se = F) + theme_bw() + scale_color_manual(values = pal) +
    xlab("Minimum intron length (bp)") + ylab("F1")
  plot.legend <- ggplot(dfpr, aes(x = min.len, y = recall, color = Tool)) + 
    geom_smooth(se = F) + theme_bw() + scale_color_manual(values = pal) +
    theme(legend.text = element_text(size = legend.text.size),
          legend.title = element_text(size = legend.text.size))
  plot.legend <- get_legend(plot.legend)
  # prep multiplot
  gg.prec <- gg.prec + theme(axis.title.x = element_blank(), legend.position = "none")
  gg.rec <- gg.rec + theme(axis.title.x = element_blank(), legend.position = "none")
  gg.f1 <- gg.f1 + theme(axis.title.x = element_blank(), legend.position = "none")
  return(list(prec = gg.prec, rec = gg.rec, f1score = gg.f1, legend = plot.legend))
})
names(lgg) <- names(ltsv)

# get formatted plot objects
prec.hx1 <- lgg[["HX1"]]$prec + xlim(0, xmax) + ggtitle("HX1") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.y = element_text(size = yaxis.title.size),
        plot.title = element_text(size = yaxis.title.size)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0, prec.ymax))

rec.hx1 <- lgg[["HX1"]]$rec + xlim(0, xmax) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.y = element_text(size = yaxis.title.size)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0, rec.ymax))

f1.hx1 <- lgg[["HX1"]]$f1score + xlim(0, xmax) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = xaxis.text.size),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.y = element_text(size = yaxis.title.size)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0, f1.ymax))


prec.ipsc <- lgg[["iPSC"]]$prec + ylim(0, prec.ymax) + 
  xlim(0, xmax) + ggtitle("iPSC") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        plot.title = element_text(size = yaxis.title.size))

rec.ipsc <- lgg[["iPSC"]]$rec + ylim(0, rec.ymax) + xlim(0, xmax) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank())

f1.ipsc <- lgg[["iPSC"]]$f1score + xlim(0, xmax) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = xaxis.text.size)) +
  ylim(0, f1.ymax)

#--------------------------------------------
# revised fig 4b -- prec/rec by ilength tiles
#--------------------------------------------
# make composite plot
# format plot vars
xlab.str <- paste0("Intron length (# bases)", paste0(rep(" ", 10), collapse = ""), collapse = "")
num.plot <- 4; num.legend <- 3
lm <- matrix(c(rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(2, num.plot + 2), rep(5, num.plot + 1), rep(7, num.legend),
               rep(2, num.plot + 2), rep(5, num.plot + 1), rep(7, num.legend),
               rep(2, num.plot + 2), rep(5, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend)),
             byrow = T, nrow = 11)

# save new plots
plot.fname <- "ggptline_prec-rec-f1_by-ilength-binned_combined-2samples"
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 7.4, 3.5)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm, 
             bottom = grid::textGrob(xlab.str,
                                     gp=grid::gpar(fontsize=title.text.size)))


dev.off()

# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 7, height = 3.5, units = "in", res = 500)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm,
             top = title.str, bottom = xlab.str)
dev.off()
