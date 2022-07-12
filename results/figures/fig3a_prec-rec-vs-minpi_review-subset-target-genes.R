#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize precision and recall, and plot with continuous features.

library(ggplot2)
library(ggpubr)
library(gridExtra)

#-----------
# load data
#-----------
# get background table
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR6026510-ipsc.csv"
tsv.fname.hx1 <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR2911306-hx1.csv"

ltsv <- list()
ltsv[["iPSC"]] <- read.csv(tsv.fname.ipsc, header = T)
ltsv[["HX1"]] <- read.csv(tsv.fname.hx1, header = T)

#-----------------
# helper functions
#-----------------
dfpr_bylrfilt <- function(tsv, lrfilt = seq(0.1, 0.9, 0.1), 
                          intron.type.cname = "filtintron", lr.metric.cname = "max_intron_persistence",
                          tool.strv = c("interest", "superintronic", "iread", "kma", "irfinders",
                                        "rmats", "majiq", "suppa2")){
  do.call(rbind, lapply(lrfilt, function(min.lr){
    do.call(rbind, lapply(tool.strv, function(tooli){
      # get main df
      tsvi <- tsv[!duplicated(tsv$intron),]
      cnv <- colnames(tsvi)
      cnv <- c(cnv[grepl(tooli, cnv)], lr.metric.cname)
      tsvi <- tsvi[,cnv]
      which.tooli.value <- grepl(".*lwm$", cnv)
      which.tooli.value <- which.tooli.value & grepl(intron.type.cname, cnv)
      tsvi[is.na(tsvi[,which.tooli.value]), which.tooli.value] <- 0
      # get prec/rec
      tp <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] > 0))
      fp <- length(which(tsvi[,lr.metric.cname] < min.lr & tsvi[,which.tooli.value] > 0))
      fn <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] == 0))
      precisioni <- tp/(tp+fp); recalli <- tp/(tp+fn)
      fm <- (2*precisioni*recalli)/(precisioni+recalli)
      data.frame(tool = tooli, min.lr = min.lr, 
                 precision = precisioni, recall = recalli, 
                 fmeasure = fm, num.lr = nrow(tsvi))
    }))
  }))
}

#----------------
# get plot data
#----------------
tool.strv = c("interest", "si", "iread", "kma", "irfinders", 
              "rmats", "majiq", "suppa2")
# get prec/rec across min lr filters
ldfpr <- lapply(ltsv, function(tsvi){
  dfpr <- dfpr_bylrfilt(tsvi)
  # get f1 score
  dfpr$`F1-score` <- (2*dfpr$precision*dfpr$recall)/(dfpr$precision+dfpr$recall) 
  # format vars
  dfpr$Tool <- ifelse(dfpr$tool=="interest", "IntEREst",
                      ifelse(dfpr$tool == "superintronic", "superintronic",
                             ifelse(dfpr$tool == "iread", "iREAD",
                                    ifelse(dfpr$tool == "kma", "KMA",
                                           ifelse(dfpr$tool == "irfinders", "IRFinder-S", 
                                                  ifelse(dfpr$tool == "rmats", "rMATS", 
                                                         ifelse(dfpr$tool == "majiq", "MAJIQ", 
                                                                ifelse(dfpr$tool == "suppa2", 
                                                                       "SUPPA2", "NA"))))))))
  dfpr
})
# names(ldfpr) <- names(ltsv)

# get the color palette:
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')
alpha.val <- 1 # set transparency
pt.size <- 1.2; line.size <- 1 # set point and line size
xlab.str <- "Minimum persistence"
lgg <- lapply(ldfpr, function(dfpr){
  # make plot objects
  gg.ptline.prec <- ggplot(dfpr, aes(x = min.lr, y = precision, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab(xlab.str) +
    ylab("Precision") + scale_color_manual(values = pal)
  gg.ptline.rec <- ggplot(dfpr, aes(x = min.lr, y = recall, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab(xlab.str) +
    ylab("Recall") + scale_color_manual(values = pal)
  gg.ptline.f1score <- ggplot(dfpr, aes(x = min.lr, y = fmeasure, color = Tool)) + 
    geom_point(alpha = alpha.val, size = pt.size) + geom_line(alpha = alpha.val, size = line.size) + 
    theme_bw() + xlab(xlab.str) +
    ylab("F1-score") + scale_color_manual(values = pal)
  plot.legend <- ggplot(dfpr, aes(x = min.lr, y = recall, color = Tool)) + 
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
names(lgg) <- c("iPSC", "HX1")

# make composite plot objects
xlim.min <- 0.1; xlim.max <- 0.9
ymax.prec <- 0.4
ymax.rec <- 0.8
ymax.f1 <- 0.28
prec.hx1 <- lgg[["HX1"]]$prec + ylim(0, ymax.prec) + xlim(xlim.min, xlim.max) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))
rec.hx1 <- lgg[["HX1"]]$rec + ylim(0, ymax.rec) + xlim(xlim.min, xlim.max) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))
f1.hx1 <- lgg[["HX1"]]$f1score + ylim(0, ymax.f1) + xlim(xlim.min, xlim.max) +
  theme(axis.title.x = element_blank()) + scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))
prec.ipsc <- lgg[["iPSC"]]$prec + ylim(0, ymax.prec) + xlim(xlim.min, xlim.max) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))
rec.ipsc <- lgg[["iPSC"]]$rec + ylim(0, ymax.rec) + xlim(xlim.min, xlim.max) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank()) +
  scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))
f1.ipsc <- lgg[["iPSC"]]$f1score + ylim(0, ymax.f1) + xlim(xlim.min, xlim.max) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = seq(xlim.min, xlim.max, 0.1))

#-----------------
# save new figures
#-----------------
num.plot <- 4; num.legend <- 2
lm <- matrix(c(rep(1, num.plot + 2), rep(4, num.plot + 1), rep(7, num.legend),
               rep(2, num.plot + 2), rep(5, num.plot + 1), rep(7, num.legend),
               rep(3, num.plot + 2), rep(6, num.plot + 1), rep(7, num.legend)),
             byrow = T, nrow = 3)
title.str <- paste0(paste0(rep(" ", 32), collapse = ""), 
                    "HX1", paste0(rep(" ", 52), collapse = ""), 
                    "iPSC", paste0(rep(" ", 90), collapse = ""),
                    collapse = "")
xlab.str <- paste0("Minimum persistence", paste0(rep(" ", 10), collapse = ""), 
                   collapse = "")

# save new figures
plot.fname <- "ggptline_prec-rec-by-lrpersistence_reviewer-subset_combined-2samples"
# pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 7, 3.5)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm,
             top = title.str, bottom = xlab.str)
dev.off()
# png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, 7, 3.5, units = "in", res = 500)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, prec.ipsc, rec.ipsc, f1.ipsc,
             as_ggplot(lgg[[1]]$legend), layout_matrix = lm,
             top = title.str, bottom = xlab.str)
dev.off()