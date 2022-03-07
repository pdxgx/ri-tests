#!/usr/bin/env R

# Author: Sean Maden
#
# Spot-check gene/RI candidates from literature review

library(ggplot2)
library(ggforce)
library(gridExtra)

#----------
# load data
#----------
titlev <- c("HX1", "iPSC")

# get background table
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "called_RI_data_summary_iPSC.tsv"
tsv.fname.hx1 <- "called_RI_data_summary_HX1.tsv"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

# get gene candidate tables
genes.fname.hx1 <- "HX1_validatedONLY_RI_genes_02-23-2022_12.56.31.tsv"
genes.fname.ipsc <- "iPSC_validatedONLY_RI_genes_02-23-2022_12.56.31.tsv"
lgenes <- list()
lgenes[["HX1"]] <- read.table(genes.fname.hx1, sep = "\t", header = T)
lgenes[["iPSC"]] <- read.table(genes.fname.ipsc, sep = "\t", header = T)

# color palette
pal <- c("FN" = "#6d9cc6", "FP" = "#d8788a", "TP" = "#9db92c", "TN" = "#f59b42")

# define params
lr.metric.cname <- "max_intron_persistence"

#-------------------------
# persistence violin plots
#-------------------------
lgg <- list()
bg.str <- "Background"
for(i in seq(length(titlev))){
  samplei <- titlev[i]
  dfp <- lgenes[[samplei]]
  dfp <- dfp[,c(lr.metric.cname, "gene_name")]
  colnames(dfp) <- c("persistence", "geneid")
  dfp.bg <- data.frame(persistence = ltsv[[samplei]][,lr.metric.cname])
  dfp.bg$geneid <- bg.str
  # order genes
  ugenev <- unique(dfp$geneid)
  ugene.medv <- unlist(lapply(ugenev, function(genei){
    median(dfp[dfp[,2]==genei & dfp[,1] > 0,1], na.rm = T)
  }))
  # order genes with bg
  dfp <- rbind(dfp.bg, dfp)
  dfp$col <- ifelse(dfp$geneid==bg.str,T,F)
  lvlv <- c(bg.str, ugenev[rev(order(ugene.medv))])
  dfp$geneid <- factor(dfp$geneid, levels = lvlv)
  # plot objects
  vp <- ggplot(dfp, aes(x = geneid, y = persistence, fill = col)) +
    geom_hline(yintercept = 0.1, color = "blue") + theme_bw() +
    geom_violin(draw_quantiles = 0.5) + ggtitle(samplei) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.position = "none")
  lgg[[samplei]] <- vp + facet_zoom(ylim = c(0, 0.15))
}

# save new figures
plot.fname <- "vp-bg-vgenes_persistence_combined-2samples"
# new pdf
pdf(paste0(plot.fname, ".pdf"), 5, 3.5)
grid.arrange(lgg[[1]], lgg[[2]], bottom = "Gene ID", left = "Persistence")
dev.off()
# new png
png(paste0(plot.fname, ".png"), width = 5, height = 3.5, units = "in", res = 500)
grid.arrange(lgg[[1]], lgg[[2]], bottom = "Gene ID", left = "Persistence")
dev.off()

#-----------------------
# truth metrics barplots
#-----------------------
# update color palette
pal <- c("FN" = "#6d9cc6", "FP" = "#d8788a", "TP" = "#9db92c")
# get plot objects
lgg <- list()
for(i in seq(length(titlev))){
  samplei <- titlev[i]; tsvi <- lgenes[[samplei]]
  tsvi <- tsvi[!duplicated(tsvi$intron),] # filter duplicate introns
  strv <- c("true_positives$", "false_negatives$", "false_positives$")
  ugenev <- unique(tsvi$gene_name)
  cnv <- colnames(tsvi)[grepl(".*true_positives$", colnames(tsvi))]
  toolv <- unique(gsub("_.*", "", cnv))
  dfp <- do.call(rbind, lapply(ugenev, function(genei){
    dfpi <- tsvi[tsvi$gene_name==genei,]
    do.call(rbind, lapply(toolv, function(tooli){
      dfpii <- dfpi[,grepl(tooli, colnames(dfpi))]
      dfprii <- do.call(rbind, lapply(strv, function(stri){
        num.value <- length(which(dfpii[,grepl(stri, colnames(dfpii))]==1))
        dfpriii <- data.frame(value = num.value)
        dfpriii$tmetric <- stri; return(dfpriii)
      }));dfprii$gene <- genei; dfprii$tool <- tooli
      return(dfprii)
    }))
  }))
  # rm nifk
  dfp <- dfp[!dfp$gene=="NIFK",]
  # format vars
  dfp$`Truth\nmetric` <- ifelse(dfp$tmetric==strv[1], "TP",
                                ifelse(dfp$tmetric==strv[2], "FN", "FP"))
  # plot data
  nrowi <- length(unique(dfp$gene))
  # counts
  bp.ct <- ggplot(dfp, aes(x = tool, y = value, fill = `Truth\nmetric`)) +
    geom_bar(stat = "identity") + ggtitle(samplei) + theme_bw() +
    scale_fill_manual(values = pal) + ylab("Introns (count)") +
    theme(axis.title.x = element_blank(), legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  bp.ct <- bp.ct + facet_wrap(~gene, nrow = nrowi)
  # percentages
  bp.perc <- ggplot(dfp, aes(x = tool, y = value, fill = `Truth\nmetric`)) +
    geom_bar(position = "fill", stat = "identity") + ggtitle("") + theme_bw() +
    scale_fill_manual(values = pal) + ylab("Introns (%)") +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme(axis.title.x = element_blank(), legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  bp.perc <- bp.perc + facet_wrap(~gene, nrow = nrowi)
  # legend object
  plot.legend <- ggplot(dfp, aes(x = tool, y = value, fill = `Truth\nmetric`)) +
    geom_bar(position = "fill", stat = "identity") + theme_bw() +
    scale_fill_manual(values = pal)
  plot.legend <- ggpubr::get_legend(plot.legend)
  # return plot list
  lgg[[samplei]][["bp.count"]] <- bp.ct
  lgg[[samplei]][["bp.perc"]] <- bp.perc
  lgg[["legend"]] <- plot.legend
}

# make new composite plot
plot.fname <- "ggbarplot-tool-tmetricct_validategenes_2samples-combined"
lm <- matrix(c(rep(c(1,2,3,4), each = 3), 5), nrow = 1, ncol = 13)
# new pdf
pdf(paste0(plot.fname, ".pdf"), 8, 6.5)
grid.arrange(lgg[["HX1"]][[1]], lgg[["HX1"]][[2]], 
             lgg[["iPSC"]][[1]], lgg[["iPSC"]][[2]],
             ggpubr::as_ggplot(lgg$legend), layout_matrix = lm)
dev.off()
# new png
png(paste0(plot.fname, ".png"), width = 8, height = 6.5, units = "in", res = 500)
grid.arrange(lgg[["HX1"]][[1]], lgg[["HX1"]][[2]], 
             lgg[["iPSC"]][[1]], lgg[["iPSC"]][[2]],
             ggpubr::as_ggplot(lgg$legend), layout_matrix = lm)
dev.off()

#-----------------------------
# truth metric calls grid plot
#-----------------------------
# update color palette
pal <- c("FN" = "#6d9cc6", "FP" = "#d8788a", "TP" = "#9db92c", "TN" = "#ffc18f")

# get the plot data
toolv <- c("IntEREst", "iREAD", "IRFinder.S", "superintronic", "kma")
cnv <- colnames(lgenes[[1]]); 
dfp <- do.call(rbind, lapply(names(lgenes), function(samplei){
  mgenei <- lgenes[[samplei]]
  do.call(rbind, lapply(toolv, function(tooli){
    cnvi <- c(cnv[grepl(tooli, cnv)], "intron", "gene_name")
    mf <- mgenei[,cnvi]; intronv <- unique(mf$intron)
    do.call(rbind, lapply(intronv, function(introni){
      mfi <- mf[mf$intron==introni,]; ni <- names(mfi)
      # message(introni)
      if(length(mfi) > 0){
        tmi <- ifelse(mfi[grepl("true_positives", ni)] == 1, "TP",
                      ifelse(mfi[grepl("false_positives", ni)] == 1, "FP", 
                             ifelse(mfi[grepl("false_negatives", ni)] == 1, "FN", "TN")))
        data.frame("tool" = tooli, "intron" = introni, "sample" = samplei,
                   "tmetric" = as.character(tmi), "gene" = mfi$gene_name)
      }
    }))
    #message(tooli)
  }))
  #message(samplei)
}))

# get plot legend
dfp$`Truth\nmetric` <- dfp$tmetric
plot.legend <- ggplot(dfp, aes(x = intron, y = tool, fill = `Truth\nmetric`)) +
  geom_tile(color = "black") + theme_bw() + scale_fill_manual(values = pal)
plot.legend <- ggpubr::get_legend(plot.legend)

# get composite plots where applicable

genev <- unique(dfp$gene)
lgg <- lapply(genev, function(genei){
  dfpi <- dfp[dfp$gene==genei,]
  samplev <- unique(dfpi$sample)
  dfpi$start <- as.numeric(gsub(".*:|-.*", "", dfpi$intron))
  dfpi$end <- as.numeric(gsub(".*-", "", dfpi$intron))
  intronv <- unique(dfpi$intron)
  startv <- as.numeric(gsub(".*:|-.*", "", intronv))
  intronv <- intronv[order(startv)]
  dfpi$intron.lab <- 0
  for(ii in seq(length(intronv))){
    dfpi[dfpi$intron==intronv[ii],]$intron.lab <- ii
  }
  
  # dfpi$intron <- factor(dfpi$intron, levels = intronv)
  ploti <- ggplot(dfpi, aes(x = intron.lab, y = tool, fill = tmetric)) +
    geom_tile(color = "black") + theme_bw() + scale_fill_manual(values = pal) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), legend.position = "none", 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
    scale_x_continuous(breaks = seq(1, max(dfpi$intron.lab), 1))
  
  ploti <- ploti + facet_wrap(~sample, nrow = length(samplev)) +
    ggtitle(genei)
  
  return(ploti)
})
names(lgg) <- genev

# save new plots
plot.fname <- "ggtile-tmetric-bytool-byintron_validate-genes_2samples-combined"
lm <- matrix(c(rep(c(1, 1, 1, 2, 2, 2, 2), 3),
               rep(c(3, 3, 3, 4, 4, 5, 5), 3),
               rep(c(6, 6, 7, 7, 8, 8, 9), 3),
               rep(11, 7)), nrow = 7, ncol = 10)
# new pdf
pdf(paste0(plot.fname, ".pdf"), 14, 9.5)
grid.arrange(lgg[["LBR"]], lgg[["CELF1"]], 
             lgg[["AP1G2"]], lgg[["IGSF8"]], lgg[["FAHD2A"]],
             lgg[["FAHD2B"]], lgg[["SRSF7"]], lgg[["CLASRP"]], lgg[["CTSD"]],
             ggpubr::as_ggplot(plot.legend),
             layout_matrix = lm, bottom = "Intron number", left = "Tool")
dev.off()
# new png
png(paste0(plot.fname, ".png"), width = 14, height = 9.5, 
    units = "in", res = 500)
grid.arrange(lgg[["LBR"]], lgg[["CELF1"]], 
             lgg[["AP1G2"]], lgg[["IGSF8"]], lgg[["FAHD2A"]],
             lgg[["FAHD2B"]], lgg[["SRSF7"]], lgg[["CLASRP"]], lgg[["CTSD"]],
             ggpubr::as_ggplot(plot.legend),
             layout_matrix = lm, bottom = "Intron number", left = "Tool")
dev.off()


