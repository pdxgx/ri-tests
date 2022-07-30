#!/usr/bin/env R

# Author: Sean Maden
#
# Analyze prec/rec results, before/after intron filt

library(ggplot2)

#----------
# load data
#----------
plot.titlev <- c("HX1", "iPSC")
tsv.fname.hx1 <- "called_RI_data_summary_HX1featureannotated_GCcontent_splicetype.tsv"
tsv.fname.ipsc <- "called_RI_data_summary_iPSCfeatureannotated_GCcontent_splicetype.tsv"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

#-----------------
# helper functions
#-----------------
# get binomial test results
do_binomtest <- function(ltsv, typev = c("HX1", "iPSC"), 
                         varnamev = c(rep("intron_type_annotation", 3), 
                                      rep("motif.binned", 3), 
                                      rep("length.binned", 2),
                                       rep("tof.binned", 2), 
                                      rep("mfb.binned", 2), 
                                      rep("bol.binned", 2)),
                         varstatv = c("u12", "u2", "other", "CTGC|GCAG", 
                                      "GTAG|CTAC", "other", "long", "short",
                                       "high", "low", "high", "low", "high", "low"),
                         tmetricv = c("tp", "fp", "fn")){
  do.call(rbind, lapply(unique(names(ltsv)), function(samplei){
    tsvi <- ltsv[[samplei]]
    lvlv <- paste0(rep(tmetricv, length(toolv)))
    lvlv <- unique(paste0(rep(varnamev, each = length(lvlv)), ";", rep(lvlv, length(varnamev))))
    dfti <- do.call(rbind, lapply(unique(lvlv), function(lvli){
      lvlvec <- unlist(strsplit(lvli,";"))
      cnvi <- colnames(tsvi)
      vari <- lvlvec[1]; tooli <- lvlvec[2]
      tmetrici <- lvlvec[3]
      which.cnvi <- which(grepl(tooli, cnvi) & grepl(tmetrici, cnvi))
      cname.testi <- cnvi[which.cnvi]
      varstatev <- unique(tsvi[,vari]) 
      varstatev <- varstatev[!is.na(varstatev)]
      # iterate on unique variable states/levels
      do.call(rbind, lapply(varstatev, function(vstatei){
        message(lvli, "; ", vstatei)
        # bg
        bg.true.cond <- tsvi[,vari] == vstatei
        bg.true.num <- length(which(bg.true.cond))
        # test
        test.total.num <- length(which(tsvi[,cname.testi] == 1))
        test.true.num <- length(which(tsvi[,vari] == vstatei & tsvi[,cname.testi] == 1))
        # binom test
        bt <- binom.test(test.true.num, test.total.num, bg.true.num/nrow(tsvi))
        bt.pval <- bt$p.value; bt.est <- bt$estimate
        data.frame(bt.pval = bt.pval, bt.est = bt.est, 
                   metric = tmetrici, variable = vari, varlevel = vstatei, 
                   sample = samplei, tool = tooli, num.bg = bg.true.num,
                   num.test = test.true.num, fract.bg = bg.true.num/nrow(tsvi),
                   fract.test = test.true.num/nrow(tsvi))
      }))
    }))
  }))
}

#-----------------------------------
# get binarized vars for stats tests
#-----------------------------------
# make the new motif var, binning identical motifs
for(i in seq(length(ltsv))){
  tsv <- ltsv[[i]]
  dt <- as.data.frame(table(tsv$motif))
  dt[,1] <- as.character(dt[,1])
  dt <- dt[rev(order(dt[,2])),]
  dt <- dt[dt[,2] > 100,] # at least 2 motif instances
  dt$revcomp <- unlist(lapply(dt[,1], function(ri){
    vect1 <- unlist(strsplit(as.character(ri), ""))
    vect2 <- ifelse(vect1 == "A", "T", 
                    ifelse(vect1 == "C", "G",
                           ifelse(vect1 == "T", "A",
                                  ifelse(vect1 == "G", "C", "NA"))))
    return(paste0(rev(vect2), collapse = ""))
  })); dt[,4] <- paste(dt[,1], dt[,3], sep = "|")
  mv <- tsv$motif
  mv <- ifelse(mv %in% c(dt[1,1], dt[1,3]), dt[1,4],
               ifelse(mv %in% c(dt[2,1], dt[2,3]), dt[2,4],
                      ifelse(mv %in% c(dt[3,1], dt[3,3]), dt[3,4],
                             ifelse(mv %in% c(dt[4,1], dt[4,3]), dt[4,4], "other"))))
  ltsv[[i]]$motif.binned <- mv
}

# intron lengths
# binned on 50th quantile
for(i in seq(length(ltsv))){
  tsvi <- ltsv[[i]]; samplei <- names(ltsv)[i]
  q50i <- quantile(tsvi$width, na.rm = T)[3]
  tsvi$length.binned <- ifelse(tsvi$width >= q50i, "long", "short")
  lbin <- list()
  for(ci.str in c("total_overlapping_features", "max_features_per_base",
                  "X._bases_overlapped", "gc_fract")){
    ci <- which(colnames(tsvi) == ci.str)
    q50i <- quantile(tsvi[tsvi[,ci] > 0,ci], na.rm = T)[3]
    lbin[[ci.str]] <- ifelse(tsvi[,ci] >= q50i, "high", "low")
  }
  tsvi$tof.binned <- lbin[[1]]; tsvi$mfb.binned <- lbin[[2]]; 
  tsvi$bol.binned <- lbin[[3]]; ltsv[[i]] <- tsvi
}

# get 4+ truth metric categories
tdfp <- do.call(rbind, lapply(names(ltsv), function(samplei){
  tsvi <- ltsv[[samplei]]; tsvi$sample <- samplei
  return(tsvi)
}))
# tp
tdfpf <- tdfp[,grepl("true_positive", colnames(tdfp))]
is.tp.four <- apply(tdfpf, 1, function(ri){length(ri[ri>0])})
tdfp$tp.four <- ifelse(is.tp.four, T, F)
table(tdfp$tp.four)
# fn
tdfpf <- tdfp[,grepl("false_negative", colnames(tdfp))]
is.fn.four <- apply(tdfpf, 1, function(ri){length(ri[ri>0])})
tdfp$fn.four <- ifelse(is.fn.four, T, F)
# fp
tdfpf <- tdfp[,grepl("false_positive", colnames(tdfp))]
is.fp.four <- apply(tdfpf, 1, function(ri){length(ri[ri>0])})
tdfp$fp.four <- ifelse(is.fp.four, T, F)

#-----------------------------------------
# binned variable barplots by sample group
#-----------------------------------------
varv <- c("motif.binned", "length.binned", "tof.binned", 
          "mfb.binned", "bol.binned", "intron_type_annotation")
dfp <- do.call(rbind, lapply(varv, function(vari){
  message(vari)
  lvlv <- unique(tdfp[,vari]); lvlv <- lvlv[!is.na(lvlv)]
  dfpi <- do.call(rbind, lapply(c("HX1", "iPSC", "background"), function(samplei){
    if(samplei=="background"){tdfpi <- tdfp[!duplicated(tdfp$intron),]
    } else{tdfpi <- tdfp[tdfp$sample==samplei,]}
    do.call(rbind, lapply(lvlv, function(lvli){
      num.tp <- nrow(tdfpi[tdfpi$tp.four == T & tdfpi[,vari]==lvli,])
      num.fn <- nrow(tdfpi[tdfpi$fn.four == T & tdfpi[,vari]==lvli,])
      num.fp <- nrow(tdfpi[tdfpi$fp.four == T & tdfpi[,vari]==lvli,])
      data.frame(tmetric = c("tp", "fn", "fp"),
                 num.vari = c(num.tp, num.fn, num.fp),
                 level = rep(lvli, 3), 
                 variable = rep(vari, 3),
                 sample = rep(samplei, 3))
    }))
  }))
}))

# format vars
dfp[dfp$sample == "background",]$sample <- "bg"
dfp$tmetric <- ifelse(dfp$tmetric=="fn", "FN", 
                      ifelse(dfp$tmetric == "fp", "FP", 
                             ifelse(dfp$tmetric == "tp", "TP", "TN")))

# get list of plots
lgg.ct <- lapply(varv, function(vari){
  bp <- ggplot(dfp[dfp$variable==vari,], 
               aes(x = sample, y = num.vari, fill = level)) + 
    geom_bar(stat = "identity") + ggtitle(paste0("Variable: ", vari)) + theme_bw() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
  bp <- bp + facet_wrap(~tmetric)
}); names(lgg.ct) <- varv
lgg.perc <- lapply(varv, function(vari){
  bp <- ggplot(dfp[dfp$variable==vari,], 
               aes(x = sample, y = 100*num.vari, fill = level)) + 
    geom_bar(stat = "identity", position = "fill") + ggtitle(" ") + theme_bw() +
    theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  bp <- bp + facet_wrap(~tmetric) + ggtitle(" ")
}); names(lgg.perc) <- varv

num.ct <- 3
num.perc <- 4
lm <- matrix(c(rep(1, num.ct), rep(2, num.perc),
               rep(3, num.ct), rep(4, num.perc),
               rep(5, num.ct), rep(6, num.perc),
               rep(7, num.ct), rep(8, num.perc),
               rep(9, num.ct), rep(10, num.perc),
               rep(11, num.ct), rep(12, num.perc)),
             byrow = T, ncol = num.ct+num.perc)

# save new plots
plot.fname <- "bp-6properties-binned-bygroup_2samples-combined"
# save new pdf
pdf(paste0(plot.fname, ".pdf"), 8, 10)
grid.arrange(lgg.ct[[1]], lgg.perc[[1]],
             lgg.ct[[2]], lgg.perc[[2]],
             lgg.ct[[3]], lgg.perc[[3]],
             lgg.ct[[4]], lgg.perc[[4]],
             lgg.ct[[5]], lgg.perc[[5]],
             lgg.ct[[6]], lgg.perc[[6]],
  layout_matrix = lm, bottom = "Sample group", left = "Number of introns")
dev.off()
# save new png
png(paste0(plot.fname, ".png"), width = 8, height = 10, units = "in", res = 500)
grid.arrange(lgg.ct[[1]], lgg.perc[[1]],
             lgg.ct[[2]], lgg.perc[[2]],
             lgg.ct[[3]], lgg.perc[[3]],
             lgg.ct[[4]], lgg.perc[[4]],
             lgg.ct[[5]], lgg.perc[[5]],
             lgg.ct[[6]], lgg.perc[[6]],
             layout_matrix = lm, bottom = "Sample group", left = "Number of introns")
dev.off()
