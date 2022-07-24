#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize precision and recall, and plot with continuous features.

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
# tsv.fname <- "nonzero_RI_data_summary_HX1featureannotated.tsv"
# tsv <- read.table(tsv.fname, sep = "\t", header = T)
# tsv.fname <- "target_genes_LR_annotated_granges-lrmap_sr-5-methods_SRR2911306-hx1.csv"
tsv.fname <- "subset_target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR2911306-hx1.csv"
tsv <- read.csv(tsv.fname, header = T)

#-----------------
# helper functions
#-----------------
# get precision and recall from a dfp of truth metrics
dfpr_bylrfilt <- function(tsv, lrfilt = seq(0.1, 0.9, 0.1), 
                          intron.type.cname = "filtintron", 
                          lr.metric.cname = "max_intron_persistence",
                          tool.strv = c("interest", "superintronic", "iread", "kma", 
                                        "irfinder", "majiq", "rmats", "suppa2")){
  do.call(rbind, lapply(lrfilt, function(min.lr){
    do.call(rbind, lapply(tool.strv, function(tooli){
      # get main df
      tsvi <- tsv[!duplicated(tsv$intron),]
      cnv <- colnames(tsvi)
      cnv <- c(cnv[grepl(tooli, tolower(cnv))], lr.metric.cname)
      tsvi <- tsvi[,cnv]
      which.tooli.value <- grepl(".*lwm$", cnv) & grepl(intron.type.cname, cnv)
      tsvi[is.na(tsvi[,which.tooli.value]), which.tooli.value] <- 0
      # get prec/rec
      tp <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] > 0))
      fp <- length(which(tsvi[,lr.metric.cname] < min.lr & tsvi[,which.tooli.value] > 0))
      fn <- length(which(tsvi[,lr.metric.cname] >= min.lr & tsvi[,which.tooli.value] == 0))
      precisioni <- tp/(tp+fp); recalli <- tp/(tp+fn)
      fm <- (2*precisioni*recalli)/(precisioni+recalli)
      tooli <- ifelse(tooli == "irfinder", "irfinders", tooli)
      data.frame(tool = tooli, min.lr = min.lr, 
                 precision = precisioni, recall = recalli, 
                 fmeasure = fm, num.lr = nrow(tsvi))
    }))
  }))
}


#-----------------------------------
# get prec/rec across min lr filters
#-----------------------------------
# get list of dfpr objects
tool.strv = c("interest", "superintronic", "iread", "kma", "irfinders", "majiq", 
              "rmats", "suppa2")
dfpr <- dfpr_bylrfilt(tsv)

# get f1 score
dfpr$`F1-score` <- (2*dfpr$precision*dfpr$recall)/(dfpr$precision+dfpr$recall)

# format vars
dfpr$Tool <- ifelse(dfpr$tool=="interest", "IntEREst",
                    ifelse(dfpr$tool == "superintronic", "superintronic",
                           ifelse(dfpr$tool == "iread", "iREAD",
                                  ifelse(dfpr$tool == "kma", "KMA",
                                         ifelse(dfpr$tool == "irfinder", "IRFinder-S", 
                                                ifelse(dfpr$tool=="majiq", "MAJIQ",
                                                       ifelse(dfpr$tool == "rmats", "rMATS",
                                                              ifelse(dfpr$tool == "suppa2", "SUPPA2", NA))))))))
#-------------------
# get quantiles data
#-------------------
# get plot data
dfp <- do.call(rbind, lapply(tool.strv, function(tooli){
  message(tooli)
  dfpri <- dfpr[dfpr$tool==tooli,]
  qi.prec <- quantile(dfpri$precision, na.rm = T)
  qi.rec <- quantile(dfpri$recall, na.rm = T)
  qi.fm <- quantile(dfpri$fmeasure, na.rm = T)
  data.frame(prec25 = qi.prec[2], prec50 = qi.prec[3], prec75 = qi.prec[4], 
             rec25 = qi.rec[2], rec50 = qi.rec[3], rec75 = qi.rec[4], 
             fm25 = qi.fm[2], fm50 = qi.fm[3], fm75 = qi.fm[4], 
             tool = tooli)
}))

# format vars
dfp[dfp$tool=="interest",]$tool <- "IntEREst"
dfp[dfp$tool=="iread",]$tool <- "iREAD"
dfp[dfp$tool=="irfinders",]$tool <- "IRFinder-S"
dfp[dfp$tool=="kma",]$tool <- "KMA"
dfp[dfp$tool=="majiq",]$tool <- "MAJIQ"
dfp[dfp$tool=="rmats",]$tool <- "rMATS"
dfp[dfp$tool=="suppa2",]$tool <- "SUPPA2"
dfp$`RI detection\ntool` <- dfp$tool

#--------------------------------------------
# composite plots -- prec, rec, f1 score iqr
#--------------------------------------------
# get the color palette:
pal <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')

# format vars
# barplots fmeasure
dfp$fm50 <- round(as.numeric(as.character(dfp$fm50)), 3)
dfp$`RI detection tool` <- dfp$`RI detection\ntool`
dfp$`RI detection tool` <- factor(dfp$`RI detection tool`, 
                                     levels = dfp$`RI detection tool`[order(dfp$fm50)])

# get plot objects
# get scatterplot object
ggpt.filt <- ggplot(dfp, aes(x = rec50, y = prec50, color = `RI detection\ntool`)) + 
  geom_point() + theme_bw() + scale_color_manual(values = pal) +
  geom_errorbar(aes(ymin = prec25, ymax = prec75), width = 0.1) + 
  geom_errorbarh(aes(xmin = rec25, xmax = rec75), height = 3e-3) + 
  xlab("Recall") + ylab("Precision") + ggtitle(plot.title)
# get barplot object
ggbp <- ggplot(dfp, aes(x = `RI detection tool`, y = fm50, fill = `RI detection\ntool`)) +
  geom_bar(stat = "identity") + theme_bw() + ylab("F1-score") + 
  geom_errorbar(aes(ymin = fm25, ymax = fm75), width = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1),
        legend.position = "none") +
  scale_fill_manual(values = pal)

# get composite plot -- pt prec/rec with bp fmeasure
ggpt.filt <- ggpt.filt + theme(legend.position = "none") + ggtitle(plot.title)
ggbp <- ggbp + ggtitle("")

# save new plots
plot.fname <- paste0("gg-pt-bp_fm-prec-rec_", srrid, "-", run.handle)
# new pdf
pdf(paste0(plot.fname, ".pdf"), 4.3, 2.2)
grid.arrange(ggpt.filt, ggbp, nrow = 1)
dev.off()
# new png
png(paste0(plot.fname, ".png"), width = 4.3, height = 2.2, units = "in", res = 500)
grid.arrange(ggpt.filt, ggbp, nrow = 1)
dev.off()
