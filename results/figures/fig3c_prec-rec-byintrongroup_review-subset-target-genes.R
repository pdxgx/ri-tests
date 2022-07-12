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
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR6026510-ipsc.csv"
tsv.fname.hx1 <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR2911306-hx1.csv"

tsv.ipsc <- read.csv(tsv.fname.ipsc)
tsv.hx1 <- read.csv(tsv.fname.hx1)
# define the lr metric
lr.metric.cname = "max_intron_persistence"
# get filtintron list
cnv.filt <- colnames(tsv.ipsc)
cnv <- c(cnv.filt[grepl("filtintron", cnv.filt)], lr.metric.cname, "intron")
ltsv <- list("iPSC" = tsv.ipsc[,cnv], "HX1" = tsv.hx1[,cnv])
# get allintron list
cnv.filt <- colnames(tsv.ipsc)
cnv <- c(cnv.filt[grepl("allintron", cnv.filt)], lr.metric.cname, "intron")
ltsv.all <- list("iPSC" = tsv.ipsc[,cnv], "HX1" = tsv.hx1[,cnv])

#-----------------
# helper functions
#-----------------
# get precision and recall from a dfp of truth metrics
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
      which.tooli.value <- grepl(".*weighted_median$|.*lwm$", cnv)
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

#-------------------
# get quantiles data
#-------------------
# get prec/rec across min lr filters
tool.strv = c("interest", "superintronic", "iread", "kma", "irfinders", "rmats",
              "majiq", "suppa2")
lltsv <- list("potential" = ltsv.all, "called" = ltsv)
ldfp <- lapply(seq(length(lltsv)), function(ii){
  intron.type <- names(lltsv)[ii]
  ltsvii <- lltsv[[ii]]
  lapply(seq(length(ltsvii)), function(jj){
    sample.type <- names(ltsvii)[jj]
    tsvjj <- ltsvii[[jj]]
    itype.cname.filt <- ifelse(intron.type == "potential", "allintron", "filtintron")
    dfpr <- dfpr_bylrfilt(tsvjj, intron.type.cname = itype.cname.filt)
    dfpr$`F1-score` <- (2*dfpr$precision*dfpr$recall)/(dfpr$precision+dfpr$recall)
    dfpr$type <- intron.type; dfpr$sample <- sample.type; return(dfpr)
  })
})
ldfp <- list("potential_ipsc" = ldfp[[1]][[2]], "potential_hx1" = ldfp[[1]][[1]],
             "called_ipsc" = ldfp[[2]][[2]], "called_hx1" = ldfp[[2]][[1]])

# get plot data
dfp <- do.call(rbind, lapply(seq(length(ldfp)), function(ii){
  dfpi <- ldfp[[ii]]
  do.call(rbind, lapply(tool.strv, function(tooli){
    dfpii <- dfpi[dfpi$tool==tooli,]
    qi.prec <- quantile(dfpii$precision); qi.rec <- quantile(dfpii$recall)
    qi.fm <- quantile(dfpii$fmeasure, na.rm = T)
    metric <- rep(c("precision", "recall", "f1score"), each = 1)
    data.frame(value50 = c(qi.prec[3], qi.rec[3], fm50 = qi.fm[3]), 
               value25 = c(qi.prec[2], qi.rec[2], qi.fm[2]), 
               value75 = c(qi.prec[4], qi.rec[4], qi.fm[4]),
               metric = metric, type = names(ldfp)[ii], tool = tooli)
  }))
}))

# format vars
dfp[dfp$tool=="interest",]$tool <- "IntEREst"
#dfp[dfp$tool=="si",]$tool <- "superintronic"
dfp[dfp$tool=="iread",]$tool <- "iREAD"
dfp[dfp$tool=="irfinders",]$tool <- "IRFinder-S"
dfp[dfp$tool=="kma",]$tool <- "KMA"
dfp[dfp$tool=="majiq",]$tool <- "MAJIQ"
dfp[dfp$tool=="rmats",]$tool <- "rMATS"
dfp[dfp$tool=="suppa2",]$tool <- "SUPPA2"
dfp$Tool <- dfp$tool
dfp$sample <- gsub(".*_", "", dfp$type)
dfp$intron_type <- gsub("_.*", "", dfp$type)

# get data by sample
# hx1
dfp.hx1 <- dfp[dfp$sample == "hx1",]
dfp.hx1.call <- dfp.hx1[dfp.hx1$intron_type == "called",]
dfp.hx1.pot <- dfp.hx1[dfp.hx1$intron_type == "potential",]
colnames(dfp.hx1.call) <- paste0(colnames(dfp.hx1.call), "_called")
colnames(dfp.hx1.pot) <- paste0(colnames(dfp.hx1.pot), "_potential")
dfp.hx1 <- cbind(dfp.hx1.pot, dfp.hx1.call)
# ipsc
dfp.ipsc <- dfp[dfp$sample == "ipsc",]
dfp.ipsc.call <- dfp.ipsc[dfp.ipsc$intron_type == "called",]
dfp.ipsc.pot <- dfp.ipsc[dfp.ipsc$intron_type == "potential",]
colnames(dfp.ipsc.call) <- paste0(colnames(dfp.ipsc.call), "_called")
colnames(dfp.ipsc.pot) <- paste0(colnames(dfp.ipsc.pot), "_potential")
dfp.ipsc <- cbind(dfp.ipsc.call, dfp.ipsc.pot)

#-----------------
# get plot objects
#-----------------
# get the color palette:
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')

# axis maxima
prec.called.max <- 0.24
prec.potential.max <- 0.13
rec.called.max <- 0.65
rec.potential.max <- 0.78
f1.called.max <- 0.25
f1.potential.max <- 0.2

# error bar widths/heights
prec.ebw.xaxis <- 0.005
prec.ebw.yaxis <- 0.0025
f1.ebw.xaxis <- 0.010 
f1.ebw.yaxis <- 0.005
rec.ebw.xaxis <- 0.04
rec.ebw.yaxis <- 0.08

prec.hx1 <- ggplot(dfp.hx1[dfp.hx1$metric_called == "precision",], 
       aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called), 
                width = prec.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = prec.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") + 
  # ylim(0, prec.called.max) + xlim(0, prec.potential.max) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank()) +
  ggtitle("Precision") +
  scale_y_continuous(limits = c(0, prec.called.max), 
                     breaks = seq(0, prec.called.max, by = 0.05)) +
  scale_x_continuous(limits = c(0, prec.potential.max), 
                     breaks = seq(0, prec.potential.max, by = 0.05))

prec.ipsc <- ggplot(dfp.ipsc[dfp.ipsc$metric_called == "precision",], 
                    aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called), 
                width = prec.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = prec.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") + 
  #ylim(0, prec.called.max) + xlim(0, prec.potential.max) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, prec.called.max), 
                     breaks = seq(0, prec.called.max, by = 0.05)) +
  scale_x_continuous(limits = c(0, prec.potential.max), 
                     breaks = seq(0, prec.potential.max, by = 0.05))

rec.hx1 <- ggplot(dfp.hx1[dfp.hx1$metric_called == "recall",], 
                  aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called), 
                width = rec.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = rec.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") + 
  ylim(0, rec.called.max) + xlim(0,rec.potential.max) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("Recall") +
  scale_y_continuous(limits = c(0, rec.called.max), 
                     breaks = seq(0, rec.called.max, by = 0.2))

rec.ipsc <- ggplot(dfp.ipsc[dfp.ipsc$metric_called == "recall",], 
                  aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called), 
                width = rec.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = rec.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") + 
  ylim(0, rec.called.max) + xlim(0,rec.potential.max) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, rec.called.max), 
                     breaks = seq(0, rec.called.max, by = 0.2))

f1.hx1 <- ggplot(dfp.hx1[dfp.hx1$metric_called == "f1score",], 
                 aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called),
                width = f1.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = f1.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") +
  # ylim(0, f1.called.max) + xlim(0, f1.potential.max) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("F1-score") +
  scale_y_continuous(limits = c(0, f1.called.max), 
                     breaks = seq(0, f1.called.max, by = 0.05)) +
  scale_x_continuous(limits = c(0, f1.potential.max), 
                     breaks = seq(0, f1.potential.max, by = 0.05))
  

f1.ipsc <- ggplot(dfp.ipsc[dfp.ipsc$metric_called == "f1score",], 
                 aes(x = value50_potential, y = value50_called, color = tool_called)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + 
  geom_errorbar(aes(ymin = value25_called, ymax = value75_called), 
                width = f1.ebw.yaxis) + 
  geom_errorbarh(aes(xmin = value25_potential, xmax = value75_potential), 
                 height = f1.ebw.xaxis) + 
  geom_abline(slope = 1, color = "black") + 
  # ylim(0, f1.called.max) + xlim(0, f1.potential.max) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, f1.called.max), 
                     breaks = seq(0, f1.called.max, by = 0.05)) +
  scale_x_continuous(limits = c(0, f1.potential.max), 
                     breaks = seq(0, f1.potential.max, by = 0.05))
  
dfp.ipsc$Tool <- dfp.ipsc$tool_called
plot.legend <- ggplot(dfp.ipsc, aes(x = value50_potential, 
                                    y = value50_called, color = Tool)) + 
  theme_bw() + geom_point() + scale_color_manual(values = pal) + theme_bw()
plot.legend <- get_legend(plot.legend)


#----------------
# save new figure
#----------------
nplot <- 3
num.legend <- 2
lm <- matrix(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(7, num.legend),
               rep(4, 3), rep(5, 3), rep(6, 3), rep(7, num.legend)), 
             nrow = 2, byrow = T)

ylab.str <- paste0(paste0(rep(" ", 5), collapse = ""), "iPSC", 
                   paste0(rep(" ", 25), collapse = ""), "HX1",
                   "\nCalled RIs", collapse = "")
xlab.str <- paste0("Potential RIs", paste0(rep(" ", 20), collapse = ""), collapse = "")

# save figure
plot.fname <- "ggpt-iqr_prec-rec-f1_reviewer-subset_combined-2samples"
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 7.2, 3.5)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, 
             prec.ipsc, rec.ipsc, f1.ipsc, 
             plot.legend, layout_matrix = lm,
             left = ylab.str, bottom = xlab.str)
dev.off()
# png fname
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 7.2, height = 3.5, units = "in", res = 500)
grid.arrange(prec.hx1, rec.hx1, f1.hx1, 
             prec.ipsc, rec.ipsc, f1.ipsc, 
             plot.legend, layout_matrix = lm,
             left = ylab.str, bottom = xlab.str)
dev.off()