#!/usr/bin/env R

# Author: Sean Maden
#
# Heatmap of prec/rec correlations with features, cutoffs, etc.

library(ggplot2)
library(ggpubr)
library(gridExtra)

srrid <- "SRR2911306"
run.handle <- "hx1"
plot.title <- "HX1"

#----------
# load data
#----------
# tsv called RIs
tsv.fname <- "called_RI_data_summary_HX1featureannotated.tsv"
tsv.bi <- read.table(tsv.fname, sep = "\t", header = T)
dim(tsv.bi)

#-----------------
# helper functions
#-----------------
# get precision and recall from a dfp of truth metrics
dfpr_byft <- function(tsv, lr.metric.cname = "max_intron_persistence",
                      tool.strv = c("interest", "superintronic", "iread", "kma", "irfinders"),
                      varv = c("width", "total_overlapping_features", "max_features_per_base", 
                               "X._bases_overlapped", "intron_position_in_tx")){
  do.call(rbind, lapply(varv, function(vari){
    message(vari)
    vari.filt <- quantile(as.numeric(tsv[tsv[,vari] > 0,vari]), 
                          seq(0,0.9,0.1), na.rm = T)
    if(vari=="width"){vari.filt <- vari.filt[1:6]}
    do.call(rbind, lapply(vari.filt, function(vari.filti){
      do.call(rbind, lapply(tool.strv, function(tooli){
        cnv <- c(colnames(tsv)[grepl(tooli, colnames(tsv))], 
                 lr.metric.cname, "width", vari)
        tsvi <- tsv[!duplicated(tsv$intron),cnv] # rm duplicate introns
        tooli.value.filt <- grepl(".*weighted_median$|.*wwm$", cnv)
        which.tooli.value <- which(tooli.value.filt)
        tsvi[is.na(tsvi[,which.tooli.value]), which.tooli.value] <- 0
        tsvi <- tsvi[tsvi[,vari] >= vari.filti,]
        which.tooli.value <- grepl(".*weighted_median$|.*lwm$", cnv)
        pos.filt <- tsvi[,lr.metric.cname] >= 0.1
        neg.filt <- tsvi[, lr.metric.cname] < 0.1
        # get truth metrics
        tp <- length(which(pos.filt & tsvi[,which.tooli.value] > 0))
        fp <- length(which(neg.filt & tsvi[,which.tooli.value] > 0))
        fn <- length(which(pos.filt & tsvi[,which.tooli.value] == 0))
        precisioni <- tp/(tp+fp)
        recalli <- tp/(tp+fn)
        fm <- (2*precisioni*recalli)/(precisioni+recalli)
        data.frame(tool = tooli, variable = vari, min.var.filt = vari.filti, 
                   precision = precisioni, recall = recalli, f1score = fm, 
                   num.lr = nrow(tsvi))
      }))
    }))
  }))
}

# get the mcor for heatmap
get_mcor <- function(tsv, round.digits = 2){
  dfft <- dfpr_byft(tsv)
  mmat <- do.call(rbind, lapply(unique(dfft$variable), function(vi){
    dffi <- dfft[dfft$variable == vi,]
    do.call(rbind, lapply(unique(dffi$tool), function(tooli){
      dffii <- dffi[dffi$tool==tooli,]
      filt.dffii <- !(dffii$precision==0 | dffii$recall==0 | dffii$f1score==0 |
                        is.na(dffii$precision) | is.na(dffii$recall) | is.na(dffii$f1score))
      dffii <- dffii[which(filt.dffii),]
      prec.ct <- cor.test(dffii$min.var.filt, dffii$precision, method = "spearman")
      rec.ct <- cor.test(dffii$min.var.filt, dffii$recall, method = "spearman")
      f1.ct <- cor.test(dffii$min.var.filt, dffii$f1score, method = "spearman")
      matrix(c(tooli, vi, prec.ct$estimate, prec.ct$p.value, rec.ct$estimate, 
               rec.ct$p.value, f1.ct$estimate, f1.ct$p.value), nrow = 1)
    }))
  }))
  mcor <- as.data.frame(mmat, stringsAsFactors = F)
  colnames(mcor) <- c("Tool", "Variable", "precision_rho", "precision_pval",
                      "recall_rho", "recall_pval", "f1score_rho", "f1score_pval")
  for(c in 3:8){mcor[,c] <- as.numeric(mcor[,c])}
  # format variables and states
  mcor[mcor$Tool=="interest",]$Tool <- "IntEREst"
  mcor[mcor$Tool=="iread",]$Tool <- "iREAD"
  mcor[mcor$Tool=="kma",]$Tool <- "KMA"
  mcor[mcor$Tool=="irfinders",]$Tool <- "IRFinder-S"
  # make the dfp for plotting
  dfp <- do.call(rbind, lapply(c("precision", "recall", "f1score"), function(mi){
    data.frame(tool = mcor$Tool, variable = mcor$Variable, metric = rep(mi, nrow(mcor)),
               value = mcor[,which(grepl(mi, colnames(mcor)))[1]])
  }))
  # format var levels
  # levels for variable
  dfp[dfp$variable == "width",]$variable <- "Intron length (bp)"
  dfp[dfp$variable == "max_features_per_base",]$variable <- "Max ft. per base"
  dfp[dfp$variable == "X._bases_overlapped",]$variable <- "Bases with overlapping ft."
  dfp[dfp$variable == "total_overlapping_features",]$variable <- "Total overlapping ft."
  # dfp[dfp$variable == "persistence",]$variable <- "Persistence"
  dfp[dfp$variable == "intron_position_in_tx",]$variable <- "Intron 5'->3' position in tx."
  # remove persistence
  dfp <- dfp[!dfp$variable == "persistence",]
  # levels for metric
  dfp[dfp$metric == "precision",]$metric <- "Precision"
  dfp[dfp$metric == "recall",]$metric <- "Recall"
  dfp[dfp$metric == "f1score",]$metric <- "F1-score"
  dfp$label <- round(dfp$value, round.digits)
  return(list(dfft = dfft, mcor = mcor, dfp = dfp))
}

#---------------------
# make/save new figure
#---------------------
lmcor <- get_mcor(tsv.bi)
dfp <- lmcor$dfp
dfp$Rho <- dfp$value

gghm <- ggplot(dfp, aes(x = tool, y = variable, fill = Rho)) + geom_tile(color = "black", size = 0.5) +
  theme_bw() + geom_text(aes(x = tool, y = variable, label = label), color = "#6B6B6B") + ylab("Intron property") +
  scale_fill_gradient2(low = "#520B71", mid = "white", high = "#E68400", name = expression(rho)) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1),
        panel.background = element_rect(fill='white', colour='black'),
        strip.background = element_rect(fill='firebrick', colour='white'),
        strip.text = element_text(color = "white"))

# make new plots
plot.fname <- paste0("hmcor_iproperty-vs-tool_",srrid,"-",run.handle)
# make new pdf
pdf(paste0(plot.fname,".pdf"), 9.2, 2.6)
gghm + facet_wrap(~metric, nrow = 1) + ggtitle(plot.title)
dev.off()
# make new png
png(paste0(plot.fname,".png"), width = 9.2, height = 2.6, units = "in", res = 500)
gghm + facet_wrap(~metric, nrow = 1) + ggtitle(plot.title)
dev.off()

