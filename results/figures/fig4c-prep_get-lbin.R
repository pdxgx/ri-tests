#!/usr/bin/env R

# Author: Sean Maden
#
# Get binned coverage summaries for intron sets

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(Rsamtools)
library(rtracklayer)
library(ggforce)
library(GenomicAlignments)

#----------
# load data
#----------
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "called_RI_data_summary_iPSC.tsv"
tsv.fname.hx1 <- "called_RI_data_summary_HX1.tsv"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

# sr bam paths
# hx1 bam paths
bam.sr.fpath.hx1 <- file.path("SRR2911306_hx1", "SRR2911306.sorted.bam")
# ipsc bam paths
bam.sr.fpath.ipsc <- file.path("SRR6026510_ipsc", "SRR6026510.sorted.bam")

#-----------------
# helper functions
#-----------------
# get intronv by group label
get_singlegroup <- function(tsv, min.tmetric = 4, regex = "true_positives$"){
  tsvf <- tsv[,c(colnames(tsv)[grepl(regex, colnames(tsv))], "intron")]
  group.data <- tsvf[apply(tsvf, 1, function(ri){length(ri[ri==1]) >= min.tmetric}),]$intron
  return(group.data)
}

# get intronv groups 
get_lgroup <- function(tsv, min.tmetric = 4,
                       regexv = c("true_positives$", "false_negatives$", "false_positives$")){
  lgroup <- list()
  grp.namev <- paste0(gsub("\\$", "", regexv), "_", min.tmetric)
  lgroup <- lapply(regexv, function(regex){
    get_singlegroup(tsv, regex, min.tmetric = min.tmetric)})
  names(lgroup) <- grp.namev
  return(lgroup)
}

# get coverages from bam file query
query_bams <- function(intronv, bam.sr.fpath){
  t1 <- Sys.time()
  ldfp <- lapply(seq(length(intronv)), function(ii){
    message("Working on intron ",ii, " of ", length(intronv), 
            ", time: ",Sys.time()-t1,"...")
    # query regions
    introni.str <- intronv[ii]
    regioni <- data.frame(chr = gsub(":.*", "", introni.str),
                          start = gsub(".*:|-.*", "", introni.str),
                          end = gsub(".*-", "", introni.str))
    gri <- makeGRangesFromDataFrame(regioni)
    # get coverages
    param <- ScanBamParam(what = c("pos", "qwidth"), which = gri)
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    # short read cov
    bamFile.sr <- bam.sr.fpath; x.sr <- scanBam(bamFile.sr, param=param)[[1]]
    cov.sr <- coverage(IRanges(x.sr[["pos"]], width=x.sr[["qwidth"]]))
    lcov <- cov.sr
    # basic coverage plots
    lengthv <- lcov@lengths
    dfpi <- do.call(rbind, lapply(seq(lengthv), function(ii){
      data.frame(start = sum(lengthv[which(seq(lengthv) < ii)]), 
                 end = sum(lengthv[which(seq(lengthv) < ii)]) + 
                   lengthv[ii])}))
    dfpi$width <- lengthv; dfpi$value <- lcov@values
    message("Finished intron ",ii, " of ", length(intronv), 
            ", time: ",Sys.time()-t1,"...")
    return(dfpi)
  }); names(ldfp) <- intronv
  return(ldfp)
}

#---------------
# get intron ids
#---------------
set.seed(0)
num.intron.pergroup <- 50
num.bins <- 1000
num.sr.ol <- 4

lgroup <- lapply(ltsv, function(tsvi){get_lgroup(tsvi, num.sr.ol)})
names(lgroup) <- names(ltsv)

#----------------
# get bam queries
#----------------
lbam <- list("iPSC" = bam.sr.fpath.ipsc, "HX1" = bam.sr.fpath.hx1)
ldfp <- lapply(seq(length(lgroup)), function(ii){
  samplei <- names(lgroup)[ii]; lgroupi <- lgroup[[ii]]
  ldfpi <- lapply(seq(length(lgroupi)), function(jj){
    tmetric <- names(lgroupi)[jj]; lgroupj <- lgroupi[[jj]]
    query_bams(intronv = lgroupi[[jj]], bam.sr.fpath = lbam[[samplei]])
  }); names(ldfpi) <- names(lgroupi)
  return(ldfpi)
}); names(ldfp) <- names(lbam)
ldfp.fname <- "ldfp_bam-sr-intron-groups_for-smooths_combined-2samples.rda"
save(ldfp, file = ldfp.fname)

# get normalized bins
interval.size <- 1/num.bins; intv <- seq(0,1,interval.size)
lnorm <- lapply(seq(length(ldfp)), function(ii){
  ldfpi <- ldfp[[ii]]; samplei <- names(ldfp)[ii]
  dfnormi <- do.call(rbind, lapply(seq(length(ldfpi)), function(jj){
    ldfpj <- ldfpi[[jj]]
    dfi <- do.call(rbind, lapply(seq(length(ldfpj)), function(ll){
      dfl <- ldfpj[[ll]]; intronl <- names(ldfpj)[ll]; dfl$intronid <- intronl
      dfl$int_max <- int_max <- as.numeric(gsub(".*-", "", intronl))
      dfl <- dfl[dfl$value > 0,] # filt 0 values
      dfl$int.start <- dfl$start/int_max
      dfl$int.end <- dfl$end/int_max; dfl
    }))
    dfi$int.start <- dfi$start/max(dfi$end)
    dfi$int.end <- dfi$end/max(dfi$end); 
    dfi$tmetric <- names(ldfpi)[jj]
    return(dfi)
  })); dfnormi$sample <- samplei
  return(dfnormi)
}); names(lnorm) <- names(ldfp)
lnorm.fname <- "lnorm_bam-sr-intron-groups_for-smooths_combined-2samples.rda"
save(lnorm, file = lnorm.fname)

# get interval summary dfs
dfint <- do.call(rbind, lapply(seq(length(lnorm)), function(ii){
  dfi <- lnorm[[ii]]; samplei <- names(lnorm)[ii]
  dfint.ii <- do.call(rbind, lapply(unique(dfi$tmetric), function(ti){
    dfii <- dfi[dfi$tmetric == ti,]
    dfint.ti <- do.call(rbind, lapply(intv, function(inti){
      which.ranges <- which(dfii$int.start >= inti & dfii$int.end < inti+interval.size)
      dfif <- dfii[which.ranges,]; dfif <- dfif[dfif$value > 0,]
      if(nrow(dfif) < 2){medi <- meani <- sdi <- vari <- 0} else{
        medi <- median(dfif$value, na.rm = T); meani <- mean(dfif$value, na.rm = T)
        sdi <- sd(dfif$value, na.rm = T); vari <- var(dfif$value, na.rm = T)}
      data.frame(bin.min = inti, bin.max = inti+interval.size,
                 median = medi, mean = meani, sd = sdi, var = vari)
    })); dfint.ti$tmetric <- ti; return(dfint.ti)
  })); dfint.ii$sample <- samplei; return(dfint.ii)
}))
dfint.fname <- "dfint_bam-sr-intron-groups_for-smooths_combined-2samples.rda"

# format vars
dfint$tmetric.label <- "NA"
dfint[grepl("true_positive.*", dfint$tmetric),]$tmetric.label <- "TP (4+)"
dfint[grepl("false_positive.*", dfint$tmetric),]$tmetric.label <- "FP (4+)"
dfint[grepl("false_negative.*", dfint$tmetric),]$tmetric.label <- "FN (4+)"
dfint$`Truth\ncategory` <- dfint$tmetric.label

# save dfint
save(dfint, file = dfint.fname)

#--------------------------
# smooth of median coverage
#--------------------------
# load dfint
dfint.fname <- "dfint_bam-sr-intron-groups_for-smooths_combined-2samples.rda"
dfint <- get(load(dfint.fname))

# get plot object
sm.median <- ggplot(dfint, aes(x = bin.min, y = median, color = `Truth\ncategory`)) + 
  geom_smooth(se = F, alpha = 0.2) + theme_bw() +
  scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank())
sm.median <- sm.median + facet_wrap(~sample, nrow = 2) +
  ylab("Median coverage")

# save figure
plot.fname <- "ggsmooth-covmedian-ibin_combined-2samples"
# save pdf
pdf.fname <- ".pdf"
pdf(paste0(plot.fname, ".pdf"), 3.5, 2.3); sm.median; dev.off()
# save png
png(paste0(plot.fname, ".png"), width = 3.5, height = 2.3, units = "in", res = 500)
sm.median; dev.off()

#--------------------------
# get variance violin plots
#--------------------------
pal <- c("FN (4+)" = "#6d9cc6", "FP (4+)" = "#d8788a", "TP (4+)" = "#9db92c")
dfint$abs.log.var <- abs(log10(dfint$var+1))

# violin plots, by sample, absolute log10 var + 1
# hx1 vp
vp.abslogvar.hx1 <- ggplot(dfint[dfint$sample=="HX1",], 
                           aes(x = tmetric.label, y = abs.log.var, fill = `Truth\ncategory`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() + ggtitle("HX1") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal)
vp.abslogvar.hx1 <- vp.abslogvar.hx1 + facet_zoom(ylim = c(0, 1))
# ipsc vp
vp.abslogvar.ipsc <- ggplot(dfint[dfint$sample=="iPSC",], 
                           aes(x = tmetric.label, y = abs.log.var, fill = `Truth\ncategory`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() + ggtitle("iPSC") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal)
vp.abslogvar.ipsc <- vp.abslogvar.ipsc + facet_zoom(ylim = c(0, 1))
# legend
plot.legend <- ggplot(dfint, aes(x = tmetric.label, y = abs.log.var, 
                          fill = `Truth\ncategory`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  scale_fill_manual(values = pal)
plot.legend <- get_legend(plot.legend)

# save new figs
plot.fname <- "ggvp-iboncov-var_combined-2samples"
ylab.str <- "|log10(var+1)|"; xlab.str <- "Truth category"
lm <- matrix(c(rep(c(rep(1,6),rep(2,6)), 4),rep(3,12)), nrow = 12)
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 5, 2.5)
grid.arrange(vp.abslogvar.hx1, vp.abslogvar.ipsc, as_ggplot(plot.legend), 
             layout_matrix = lm, left = ylab.str)
dev.off()
# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 5, height = 2.5, res = 500, units = "in")
grid.arrange(vp.abslogvar.hx1, vp.abslogvar.ipsc, as_ggplot(plot.legend), 
             layout_matrix = lm, left = ylab.str)
dev.off()

#--------------------------------------
# ttests of variances by metric, sample
#--------------------------------------
samplev <- c("iPSC", "HX1")
for(samplei in samplev){
  dfi <- dfint[dfint$sample==samplei,]
  group1 <- c("true_positives_4", "true_positives_4", "false_negatives_4")
  group2 <- c("false_negatives_4", "false_positives_4", "false_positives_4")
  for(ii in seq(3)){
    tti <- t.test(dfi[dfi$tmetric == group1[ii],]$abs.log.var,
                  dfi[dfi$tmetric == group2[ii],]$abs.log.var)
    message("Sample: ", samplei)
    message("T-test: group1: ",group1[ii], ", group2: ", group2[ii])
    message("p-val: ", tti$p.value)
  }
}

#---------------
# parse gtf data
#---------------
# get gencode v35 gtf
# gtf source: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
gtf.fname <- "gencode.v35.annotation.gtf"
gtf <- rtracklayer::import(gtf.fname, format = "gtf")
gtf.exons <- gtf[which(gtf@elementMetadata$type=="exon")] # get exons only
# parse bins
num.bins <- 1000
lgtfi <- lapply(seq(length(ldfp)), function(samplei){
  ldfpi <- ldfp[[samplei]]; namev <- names(ldfpi)
  dftmi <- do.call(rbind, lapply(seq(length(namev)), function(namei){
    intronv <- names(ldfpi[[namei]]); tmetrici <- namev[namei]
    dfintv <- data.frame(start = gsub(".*:|-.*", "", intronv),
                         end = gsub(".*-", "", intronv),
                         chr = gsub(":.*", "", intronv))
    grintv <- makeGRangesFromDataFrame(dfintv)
    gtf.exons.int <- subsetByOverlaps(gtf.exons, grintv)
    grintvf <- subsetByOverlaps(grintv, gtf.exons.int)
    dex <- do.call(rbind, lapply(seq(length(grintvf)), function(ii){
      grintii <- grintvf[ii]
      # get exon counts
      se.ol <- summarizeOverlaps(gtf.exons.int, grintii)
      countv <- assays(se.ol)$count
      exonv <- rowRanges(se.ol)[which(countv>0)]
      if(length(exonv) > 0){
        # get intron bins
        coordv <- c(start(grintii), end(grintii))
        max.pos <- max(coordv)
        min.pos <- min.count <- min(coordv)
        bin.width <- abs(max.pos-min.pos)/num.bins
        seqv <- seq(min.pos, max.pos, bin.width)
        bins.df <- data.frame(start = seqv[1:length(seqv)-1],
                              end = seqv[2:length(seqv)])
        bins.df$chr <- unique(seqnames(grintii))
        rangesv <- gsub("\\..*", "", paste0(bins.df$start,";",bins.df$end))
        bins.df <- bins.df[!duplicated(rangesv),]; dim(bins.df)
        bins.gr <- makeGRangesFromDataFrame(bins.df)
        bins.gr@elementMetadata$bin.index <- seq(length(bins.gr))
        bins.gr@elementMetadata$exon.count <- 0
        exon.countv <- rep(0, length(bins.gr))
        for(jj in seq(length(exonv))){
          fol.gr <- findOverlaps(bins.gr, exonv[jj])
          exon.countv[queryHits(fol.gr)] <- exon.countv[queryHits(fol.gr)] + 1
        }; message(ii, "/", length(grintvf))
        dfi <- data.frame(index = seq(length(bins.gr)), exon.count = exon.countv)
        dfi <- dfi[dfi[,2] > 0,]; return(dfi)
      }
    }))
    dfex <- do.call(rbind, lapply(unique(dex[,1]), function(ii){
      valv <- dex[dex[,1]==ii,2]
      data.frame(ii, length(valv), median(valv), mean(valv), var(valv), max(valv))
    }));dfex$tmetric <- tmetrici;dfex$fract.intron <- dfex[,2]/length(intronv)
    return(dfex)}))
  colnames(dftmi) <- c("index", "num.introns.exonoverlap", "median.exons", "mean.exons", 
                       "var.exons", "max.exons", "tmetric", "fract.introns.exonoverlap")
  dftmi$sample <- names(ldfp)[samplei];return(dftmi)
})

# get plot data
dfp.gtf <- do.call(rbind, lgtfi)
# format vars
dfp.gtf$tmetric.label <- "NA"
dfp.gtf[grepl("true_positive.*", dfp.gtf$tmetric),]$tmetric.label <- "TP (4+)"
dfp.gtf[grepl("false_positive.*", dfp.gtf$tmetric),]$tmetric.label <- "FP (4+)"
dfp.gtf[grepl("false_negative.*", dfp.gtf$tmetric),]$tmetric.label <- "FN (4+)"
dfp.gtf$`Truth\ncategory` <- dfp.gtf$tmetric.label

# save data
dfp.gtf.fname <- "dfp-gtf-exonol_combined-2samples.rda"
save(dfp.gtf, file = dfp.gtf.fname)

#-------------------
# exon overlap plots
#-------------------
dfp.gtf <- get(load("dfp-gtf-exonol_combined-2samples.rda"))

# plot results -- num introns w/exon overlapping
# get plot object
gtf.smooth <- ggplot(dfp.gtf, aes(x = index, y = num.introns.exonoverlap, 
                                  color = `Truth\ncategory`)) +
  geom_smooth(se = F) + theme_bw() + ylab("Exon overlaps\n(number of introns)") +
  scale_x_continuous(breaks = c(0, 1000), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank())
gtf.smooth <- gtf.smooth + facet_wrap(~sample, nrow = 2)
# save new figure
pdf.fname <- "ggsmooth-gtf_exonoverlap-numint_combined-2samples.pdf"
pdf(pdf.fname, 3.5, 2.3); gtf.smooth; dev.off()

# plot results -- fraction introns w/exon overlapping
# get plot object
gtf.smooth <- ggplot(dfp.gtf, aes(x = index, y = fract.introns.exonoverlap, 
                                  color = `Truth\ncategory`)) +
  geom_smooth(se = F) + theme_bw() + ylab("Exon overlaps\n(fraction of introns)") +
  scale_x_continuous(breaks = c(0, 1000), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank())
gtf.smooth <- gtf.smooth + facet_wrap(~sample, nrow = 2)
# save new figure
pdf.fname <- "ggsmooth-gtf_exonoverlap-fractint_combined-2samples.pdf"
pdf(pdf.fname, 3.5, 2.3); gtf.smooth; dev.off()

# plot results -- fraction introns w/exon overlapping, rectangles
# get plot object
gtf.smooth <- ggplot(dfp.gtf, aes(xmin = index, xmax = index, ymin=0, ymax = fract.introns.exonoverlap, 
                                  color = `Truth\ncategory`)) +
  geom_rect() + theme_bw() + ylab("Exon overlaps\n(fraction of introns)") +
  scale_x_continuous(breaks = c(0, 1000), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank())
gtf.smooth <- gtf.smooth + facet_wrap(~sample, nrow = 2)
# save new figure
pdf.fname <- "ggsmooth-gtf_exonoverlap-fractint_combined-2samples.pdf"
pdf("newplot.pdf", 3.5, 2.3); 
plot <- ggplot(dfp.gtf, aes(xmin = index, xmax = index, ymin=0, ymax = fract.introns.exonoverlap, 
                            color = `Truth\ncategory`)) +
  geom_rect(alpha = 0.4)
plot + facet_wrap(~`Truth\ncategory`)
dev.off()

#------------------------
# final composite smooths
#------------------------
library(scales)
# load data
# gtf data
gtf.fname <- "dfp-gtf-exonol_combined-2samples.rda"
dfp.gtf <- get(load(gtf.fname))
# interval coverages
dfint.fname <- "dfint_bam-sr-intron-groups_for-smooths_combined-2samples.rda"
dfint <- get(load(dfint.fname))

# color palette
pal <- c("FN (4+)" = "#6d9cc6", "FP (4+)" = "#d8788a", "TP (4+)" = "#9db92c")
dfint$log.median <- log10(dfint$median + 1)
ylab.str.coverage <- "Coverage\n(log10 scale)"
ylab.str.exon <- "Fraction of introns\noverlapping exons"
# get plot objects
dfii <- dfint[dfint$sample=="HX1",]; dfii$median <- dfii$median + 1
dfii$log.median <- log10(dfii$median)
median.hx1 <- ggplot(dfii, aes(x = bin.min, y = log.median, color = `Truth\ncategory`)) + 
  geom_smooth(se = F, alpha = 0.2) + theme_bw() + ylab(paste0(ylab.str.coverage)) +
  scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.text.x = element_blank(),
        plot.margin = unit(c(1.8, 1.8, 1.8, 3), "mm")) + ggtitle("HX1") +
  scale_color_manual(values = pal)
dfii <- dfint[dfint$sample=="iPSC",]; dfii$median <- dfii$median + 1
dfii$log.median <- log10(dfii$median)
median.ipsc <- ggplot(dfii, aes(x = bin.min, y = log.median, color = `Truth\ncategory`)) + 
  geom_smooth(se = F, alpha = 0.2) + theme_bw() +
  scale_x_continuous(breaks = c(0, 1), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none", axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ggtitle("iPSC") +
  scale_color_manual(values = pal)
# gtf smooths
gtf.smooth.hx1 <- ggplot(dfp.gtf[dfp.gtf$sample == "HX1",], 
                         aes(x = index, y = fract.introns.exonoverlap, 
                             color = `Truth\ncategory`)) +
  geom_smooth(se = F) + theme_bw() + ylab(ylab.str.exon) +
  scale_x_continuous(breaks = c(0, 1000), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = pal)
gtf.smooth.ipsc <- ggplot(dfp.gtf[dfp.gtf$sample == "iPSC",], 
                          aes(x = index, y = fract.introns.exonoverlap, 
                              color = `Truth\ncategory`)) +
  geom_smooth(se = F) + theme_bw() + 
  scale_x_continuous(breaks = c(0, 1000), labels=c("5'", "3'")) +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  scale_color_manual(values = pal)

# plot legend
plot.legend <- ggplot(dfint, aes(x = bin.min, y = median, color = `Truth\ncategory`)) + 
  geom_smooth(se = F, alpha = 0.2) + theme_bw() + scale_color_manual(values = pal)
plot.legend <- get_legend(plot.legend)

# save new figure
plot.fname <- "ggsmooth-mediancov-exonolfract_combined-2samples"
lm <- matrix(c(rep(c(1,2), 10), rep(c(3,4), 8), rep(5, 6)), nrow = 2)
# save new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 5.5, 3)
grid.arrange(median.hx1, gtf.smooth.hx1, median.ipsc, gtf.smooth.ipsc,
             as_ggplot(plot.legend), layout_matrix = lm, nrow = 2)
dev.off()
# save new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 5.5, height = 3, units = "in", res = 500)
grid.arrange(median.hx1, gtf.smooth.hx1, median.ipsc, gtf.smooth.ipsc,
             as_ggplot(plot.legend), layout_matrix = lm, nrow = 2)
dev.off()

#-----------
# corr tests
#-----------
samplev <- c("HX1", "iPSC")
for(samplei in samplev){
  dfexi <- dfp.gtf[dfp.gtf$sample==samplei, 
                   c("index", "fract.introns.exonoverlap", "tmetric")]
  dfinti <- dfint[dfint$sample==samplei, 
                  c("bin.min", "log.median", "tmetric")]
  for(tmetrici in c("true_positives_4", "false_positives_4", 
                    "false_negatives_4")){
    dfintii <- dfinti[dfinti$tmetric==tmetrici,]
    dfexii <- do.call(rbind, lapply(c(0, seq(1000)), function(ii){
      dfi <- data.frame(pos = ii/1000); dfi$fi.ex <- 0
      if(ii %in% dfexi$index){
        dfi$fi.ex <- median(dfexi[dfexi$index==ii,]$fract.introns.exonoverlap)
      }
      return(dfi)
    }))
    dfexii$bin.min <- as.numeric(dfexii$pos)
    dfexii <- dfexii[order(dfexii$bin.min, dfintii$bin.min),]
    cond <- identical(as.character(dfexii$bin.min), 
                      as.character(dfintii$bin.min))
    if(cond){
      cti <- cor.test(dfexii$fi.ex, dfintii$log.median, 
                      method = "spearman")
      message("sample: ", samplei)
      message("tmetric: ", tmetrici)
      message("rho: ", round(cti$estimate, digits = 3))
      message("pval: ", format(cti$p.value, scientific = T))
    }
  }
}

# sample: HX1
# tmetric: true_positives_4
# rho: 0.198
# pval: 2.455e-10
# sample: HX1
# tmetric: false_positives_4
# rho: 0.284
# pval: 5.091e-20
# sample: HX1
# tmetric: false_negatives_4
# rho: 0.389
# pval: 1.501e-37
# sample: iPSC
# tmetric: true_positives_4
# rho: 0.12
# pval: 1.441e-04
# sample: iPSC
# tmetric: false_positives_4
# rho: 0.209
# pval: 2.372e-11
# sample: iPSC
# tmetric: false_negatives_4
# rho: 0.328
# pval: 1.632e-26
