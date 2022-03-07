#!/usr/bin/env R

# Author: Sean Maden
#
# Make pileup plots (SR and LR results) for targeted retained/non-retained introns.

library(ggplot2)
library(Rsamtools)
library(ggforce)
library(ggpubr)
library(gridExtra)

#-----------------
# helper functions
#-----------------
get_dfp <- function(lbin, save.fname){
  dfp.all <- do.call(rbind,lapply(seq(length(lbin)), function(lii){
    dfpi <- lbin[[lii]]$dfnorm; dfpi$tmetric <- names(lbin)[lii]; dfpi}))
  message("Formatting variables and variable states...")
  dfp.all$Type <- dfp.all$type
  dfp.all[dfp.all$Type == "lr",]$Type <- "LR"; dfp.all[dfp.all$Type == "sr",]$Type <- "SR"
  strv <- c("true_positives", "false_negatives", "false_positives")
  strv.new <- c("True positives", "False negatives", "False positives")
  for(ii in seq(length(strv))){
    new.str <- paste0(strv.new[ii]," (", lbin[[ii]]$num.introns," introns)")
    dfp.all[grepl(strv[ii], dfp.all$tmetric),]$tmetric <- new.str
  }
  message("Saving new data...");save(dfp.all, file = save.fname);return(dfp.all)
}

# plot params
xmin = 0
xmax = 1
aval = 0.4
lsize = 1
plot.width = 3
plot.height = 3

#--------------
# get bar plots
#--------------
srrid <- "SRR2911306"; run.handle <- "hx1"; plot.title <- "HX1"
# get plot data
lbin <- get(load(paste0('lbin-normcov_',srrid, "-", run.handle, ".rda")))
# get all plot data
dfp.fname <- paste0("dfp-lr-sr-coverage-normloc_", srrid, "-", run.handle, ".rda")
dfp.hx1 <- get_dfp(lbin, dfp.fname)
# plot objects
dfp <- dfp.hx1; dfp <- dfp[dfp$Type=="SR",]
strv <- c("True positives", "False positives", "False negatives")
lgg.hx1 <- lapply(strv, function(stri){
  which.row <- which(grepl(stri, dfp$tmetric))
  dfpi <- dfp[which.row,]; type.str <- unique(dfpi$tmetric)
  # get max values
  text.size = 2.8
  max.sr <- max(dfpi[dfpi$Type == "SR",]$median)
  offset.pos.sr <- -0.8
  # get plot object
  ggsmooth <- ggplot(data = dfpi, aes(x = bin.min, y = median)) + 
    geom_hline(yintercept = 0, color = "black") + theme_bw() + 
    geom_bar(stat = "identity", alpha = aval) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
    geom_hline(yintercept = max.sr, color = "black", alpha = 0.4,
               linetype = "dotted") +
    geom_text(x = 0.2, y = max.sr + offset.pos.sr, label = round(max.sr,2), 
              color = "black", size = text.size, alpha = 0.4)
  ggsmooth <- ggsmooth + facet_zoom(ylim = c(-1, 2)) + ggtitle(type.str)
  ggsmooth
}); names(lgg.hx1) <- strv

# get ipsc plots
srrid <- "SRR6026510"; run.handle <- "ipsc"; plot.title <- "iPSC"
# get plot data
lbin <- get(load(paste0('lbin-normcov_',srrid, "-", run.handle, ".rda")))
# get all plot data
dfp.fname <- paste0("dfp-lr-sr-coverage-normloc_", srrid, "-", run.handle, ".rda")
dfp.ipsc <- get_dfp(lbin, dfp.fname)
# plot objects
dfp <- dfp.ipsc; dfp <- dfp[dfp$Type=="SR",]
strv <- c("True positives", "False positives", "False negatives")
lgg.ipsc <- lapply(strv, function(stri){
  which.row <- which(grepl(stri, dfp$tmetric))
  dfpi <- dfp[which.row,];type.str <- unique(dfpi$tmetric)
  # get max values
  text.size = 2.8; offset.pos.sr <- -0.8; max.sr <- max(dfpi$median)
  # get plot object
  ggsmooth <- ggplot(data = dfpi, aes(x = bin.min, y = median)) + 
    geom_hline(yintercept = 0, color = "black") + theme_bw() + 
    geom_bar(stat = "identity", alpha = aval) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
    geom_hline(yintercept = max.sr, color = "black", alpha = 0.4,
               linetype = "dotted") +
    geom_text(x = 0.2, y = max.sr + offset.pos.sr, label = round(max.sr,2), 
              color = "black", size = text.size, alpha = 0.4)
  ggsmooth <- ggsmooth + facet_zoom(ylim = c(-1, 2)) + ggtitle(type.str)
  ggsmooth
}); names(lgg.ipsc) <- strv

# get plot format
num.plot <- 8; num.legend <- 2
lm <- matrix(seq(6), nrow = 3)
# save new figure
title.str <- paste0(paste0(rep(" ", 10), collapse = ""),"HX1",
                    paste0(rep(" ", 85), collapse = ""), "iPSC",
                    paste0(rep(" ", 105), collapse = ""), collapse = "")
ylab.str <- paste0(paste0(rep(" ", 1), collapse = ""), 
                   "Median binned coverage")
pdf("gg-barplot_icov-lr-sr_combined-2samples.pdf", 10, 5)
grid.arrange(lgg.hx1[[1]], lgg.hx1[[2]], lgg.hx1[[3]], 
             lgg.ipsc[[1]], lgg.ipsc[[2]], lgg.ipsc[[3]], 
             layout_matrix = lm,
             top = title.str, left = ylab.str)
dev.off()

#-------------------
# get smoothed plots
#-------------------
srrid <- "SRR2911306"; run.handle <- "hx1"; plot.title <- "HX1"
# get plot data
lbin <- get(load(paste0('lbin-normcov_',srrid, "-", run.handle, ".rda")))
# get all plot data
dfp.fname <- paste0("dfp-lr-sr-coverage-normloc_", srrid, "-", run.handle, ".rda")
dfp.hx1 <- get_dfp(lbin, dfp.fname)
# plot objects
dfp <- dfp.hx1
strv <- c("True positives", "False positives", "False negatives")
lgg.hx1 <- lapply(strv, function(stri){
  which.row <- which(grepl(stri, dfp$tmetric))
  dfpi <- dfp[which.row,]; type.str <- unique(dfpi$tmetric)
  # get max values
  text.size = 2.8
  loess.lr <- loess(median ~ bin.min, data = dfpi[dfpi$Type=="LR",])
  loess.sr <- loess(median ~ bin.min, data = dfpi[dfpi$Type=="SR",])
  max.lr <- summary(predict(loess.lr))[6]
  max.sr <- summary(predict(loess.sr))[6]
  offset.pos.inc <- 1.2
  offset.pos.lr <- ifelse(max.lr > max.sr, -1*offset.pos.inc, offset.pos.inc)
  offset.pos.sr <- ifelse(max.sr > max.lr, -1*offset.pos.inc, offset.pos.inc)
  # get plot object
  ggsmooth <- ggplot(data = dfpi, aes(x = bin.min, y = median, color = Type)) + 
    geom_hline(yintercept = 0, color = "black") + theme_bw() + 
    stat_smooth(method="loess", se = F, alpha = aval, size = lsize) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
    geom_hline(yintercept = max.lr, color = "#F8766D", alpha = 0.4) +
    geom_hline(yintercept = max.sr, color = "#00BFC4", alpha = 0.4) +
    geom_text(x = 0.2 + 0.2, y = max.lr + offset.pos.lr, label = round(max.lr,2), 
              color = "#F8766D", size = text.size, alpha = 0.4) +
    geom_text(x = 0.2, y = max.sr + offset.pos.sr, label = round(max.sr,2), 
              color = "#00BFC4", size = text.size, alpha = 0.4)
  ggsmooth <- ggsmooth + facet_zoom(ylim = c(-1, 2)) + ggtitle(type.str)
  ggsmooth
}); names(lgg.hx1) <- strv

# get ipsc plots
srrid <- "SRR6026510"; run.handle <- "ipsc"; plot.title <- "iPSC"
# get plot data
lbin <- get(load(paste0('lbin-normcov_',srrid, "-", run.handle, ".rda")))
# get all plot data
dfp.fname <- paste0("dfp-lr-sr-coverage-normloc_", srrid, "-", run.handle, ".rda")
dfp.ipsc <- get_dfp(lbin, dfp.fname)
# plot objects
dfp <- dfp.ipsc
strv <- c("True positives", "False positives", "False negatives")
lgg.ipsc <- lapply(strv, function(stri){
  which.row <- which(grepl(stri, dfp$tmetric))
  dfpi <- dfp[which.row,]; type.str <- unique(dfpi$tmetric)
  # get max values
  text.size = 2.8
  loess.lr <- loess(median ~ bin.min, data = dfpi[dfpi$Type=="LR",])
  loess.sr <- loess(median ~ bin.min, data = dfpi[dfpi$Type=="SR",])
  max.lr <- summary(predict(loess.lr))[6]
  max.sr <- summary(predict(loess.sr))[6]
  offset.pos.inc <- 0.4
  offset.pos.lr <- ifelse(max.lr > max.sr, -1*offset.pos.inc, offset.pos.inc)
  offset.pos.sr <- ifelse(max.sr > max.lr, -1*offset.pos.inc, offset.pos.inc)
  # get plot object
  ggsmooth <- ggplot(data = dfpi, aes(x = bin.min, y = median, color = Type)) + 
    geom_hline(yintercept = 0, color = "black") + theme_bw() + 
    stat_smooth(method="loess", se = F, alpha = aval, size = lsize) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'"))+
    geom_hline(yintercept = max.lr, color = "#F8766D", alpha = 0.4) +
    geom_hline(yintercept = max.sr, color = "#00BFC4", alpha = 0.4) +
    geom_text(x = 0.2 + 0.2, y = max.lr + offset.pos.lr, label = round(max.lr,2), 
              color = "#F8766D", size = text.size, alpha = 0.4) +
    geom_text(x = 0.2, y = max.sr + offset.pos.sr, label = round(max.sr,2), 
              color = "#00BFC4", size = text.size, alpha = 0.4)
   ggsmooth <- ggsmooth + facet_zoom(ylim = c(-1, 2)) + ggtitle(type.str)
  ggsmooth
}); names(lgg.ipsc) <- strv

# save new composite figure -- ggsmooth
# get plot legend
plot.legend <- ggplot(data = dfp.all, aes(x = bin.min, y = median, color = Type)) +
  stat_smooth(geom="line", se = F, alpha = aval) + theme_bw()
plot.legend <- get_legend(plot.legend)

# get plot format
num.plot <- 8; num.legend <- 2
lm <- matrix(c(rep(1, num.plot), rep(4, num.plot), rep(7, num.legend),
               rep(2, num.plot), rep(5, num.plot), rep(7, num.legend),
               rep(3, num.plot), rep(6, num.plot), rep(7, num.legend)), byrow = T, nrow = 3)
# save new figure
title.str <- paste0(paste0(rep(" ", 10), collapse = ""),"HX1",
                    paste0(rep(" ", 85), collapse = ""), "iPSC",
                    paste0(rep(" ", 105), collapse = ""), collapse = "")
ylab.str <- paste0(paste0(rep(" ", 1), collapse = ""), 
                   "Median binned coverage")
pdf("gg-smooth_icov-lr-sr_combined-2samples.pdf", 10, 5)
grid.arrange(lgg.hx1[[1]], lgg.hx1[[2]], lgg.hx1[[3]], 
             lgg.ipsc[[1]], lgg.ipsc[[2]], lgg.ipsc[[3]], 
             as_ggplot(plot.legend), layout_matrix = lm,
             top = title.str, left = ylab.str)
dev.off()

#------------------------------------------
# save new composite figure -- violin plots
#------------------------------------------
ldfp <- list(iPSC = dfp.ipsc, HX1 = dfp.hx1)
ldfp <- lapply(seq(length(ldfp)), function(ii){
  dfp <- ldfp[[ii]]
  sample.name <- names(ldfp)[ii]
  vp.col <- ifelse(sample.name == "iPSC", "navyblue", "firebrick")
  # format vars
  dfp[grepl("True ", dfp$tmetric),]$tmetric <- "TP"
  dfp[grepl("False pos.*", dfp$tmetric),]$tmetric <- "FP"
  dfp[grepl("False neg.*", dfp$tmetric),]$tmetric <- "FN"
  dfp$sample <- names(ldfp)[ii]; dfp
})

# get plot objects
# ipsc
vp.ipsc <- ggplot(ldfp[[1]], aes(x = tmetric, y = median, fill = Type)) + 
  geom_violin(draw_quantiles = 0.5, color = "black")
vp.ipsc <- vp.ipsc + facet_zoom(ylim = c(0, 5)) + theme_bw() +
  ggtitle("iPSC") + theme(axis.title.y = element_blank(), 
                          axis.title.x = element_blank(),
                          legend.position = 'none')
# hx1
vp.hx1 <- ggplot(ldfp[[2]], aes(x = tmetric, y = median, fill = Type)) + 
  geom_violin(draw_quantiles = 0.5, color = "black")
vp.hx1 <- vp.hx1 + facet_zoom(ylim = c(0, 5)) + theme_bw() +
  ggtitle("HX1") + theme(axis.title.y = element_blank(), 
                         axis.title.x = element_blank(),
                         legend.position = 'black')
# plot legend
plot.legend <- ggplot(ldfp[[2]], aes(x = tmetric, y = median, fill = Type)) + 
  geom_violin(draw_quantiles = 0.5, color = "black") + theme_bw()
plot.legend <- get_legend(plot.legend)

# save new plot
nplot <- 5
lm <- matrix(c(rep(1, nplot), 3, rep(2, nplot), 3), nrow = 2, byrow = T)
pdf.fname <- "ggvp-medianbincov-withzoom_combined-2samples.pdf"
pdf(pdf.fname, 6, 3)
grid.arrange(vp.hx1, vp.ipsc, as_ggplot(plot.legend), 
             nrow = 2, left = "Median binned coverage",
             layout_matrix = lm)
dev.off()

#---------------------------
# summary stats, stats tests
#---------------------------
typev <- c("SR", "LR"); sampv <- c("ipsc", "hx1")
for(i in seq(2)){
  dfpi <- ldfp[[i]]
  for(typei in typev){
    message(sampv[i], ": ", typei, ": TP var = ", 
            var(dfpi[dfpi$tmetric=="TP" & dfpi$Type==typei,]$median))
    message(sampv[i], ": ", typei, ": FN var = ", 
            var(dfpi[dfpi$tmetric=="FN" & dfpi$Type==typei,]$median))
  }
}
# ipsc: SR: TP var = 12.3690628128128
# ipsc: SR: FN var = 22.0158636136136
# ipsc: LR: TP var = 25.3667417417417
# ipsc: LR: FN var = 333.655511511512
# hx1: SR: TP var = 2.24336236236236
# hx1: SR: FN var = 13.4664061561562
# hx1: LR: TP var = 120.803578578579
# hx1: LR: FN var = 4175.30968468468

# test differential variances
dfpi <- ldfp[[1]]
vt1 <- var.test(dfpi[dfpi$tmetric=="TP" & dfpi$Type=="SR",]$median, 
         dfpi[dfpi$tmetric=="FN" & dfpi$Type=="SR",]$median)
vt2 <- var.test(dfpi[dfpi$tmetric=="TP" & dfpi$Type=="LR",]$median, 
         dfpi[dfpi$tmetric=="FN" & dfpi$Type=="LR",]$median)

dfpi <- ldfp[[2]]
vt3 <- var.test(dfpi[dfpi$tmetric=="TP" & dfpi$Type=="SR",]$median, 
         dfpi[dfpi$tmetric=="FN" & dfpi$Type=="SR",]$median)
vt4 <- var.test(dfpi[dfpi$tmetric=="TP" & dfpi$Type=="LR",]$median, 
         dfpi[dfpi$tmetric=="FN" & dfpi$Type=="LR",]$median)

min(c(vt1$estimate, vt2$estimate, vt3$estimate, vt4$estimate)) # 0.02893284
max(c(vt1$estimate, vt2$estimate, vt3$estimate, vt4$estimate)) # 0.561825
summary(c(vt1$p.value, vt2$p.value, vt3$p.value, vt4$p.value)) 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 0.000e+00 0.000e+00 3.656e-20 3.656e-20 1.462e-19
