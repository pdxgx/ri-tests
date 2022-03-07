#!/usr/bin/env R

# Author: Sean Maden
# 
# Analyze the short read data, mapped to long read intron ranges.

library(ggplot2)
library(UpSetR)

srrid <- "SRR6026510"; run.handle <- "ipsc"

# load data
irv.fname <- paste0("granges-lrmap_sr-4-methods_",
                    srrid,"-",run.handle,".rda")
irv <- get(load(irv.fname))

#-------------
# correlations
#-------------
library(ggplot2)
library(tidyverse)

# get metadata df
dati <- as.data.frame(irv@elementMetadata, stringsAsFactors = F)
which.na <- is.na(dati$kma_retention_weighted_median)
dati[which.na,]$kma_retention_weighted_median <- 0

# get correlation matrix
mcor <- cor(dati[,grepl("_weighted_median$", colnames(dati))], method = "spearman")

# get the dfp plot data
dfp <- mcor
which.na <- which(upper.tri(dfp)|dfp==1); dfp[which.na] <- "NA"
dfp <- reshape2::melt(dfp)
dfp$value <- as.numeric(dfp$value)
dfp$label <- as.character(round(dfp$value, digits = 2))
dfp$Var1 <- as.character(dfp$Var1)
dfp$Var2 <- as.character(dfp$Var2)

# format variable names
dfp[dfp$Var1 == "iread_fpkm_weighted_median",]$Var1 <- "iREAD (FPKM)"
dfp[dfp$Var2 == "iread_fpkm_weighted_median",]$Var2 <- "iREAD (FPKM)"
dfp[dfp$Var1 == "interest_intret_freq_weighted_median",]$Var1 <- "IntEREst (IntRet)"
dfp[dfp$Var2 == "interest_intret_freq_weighted_median",]$Var2 <- "IntEREst (IntRet)"
dfp[dfp$Var1 == "si_score_weighted_median",]$Var1 <- "superintronic (score)"
dfp[dfp$Var2 == "si_score_weighted_median",]$Var2 <- "superintronic (score)"
dfp[dfp$Var1 == "kma_tpm_weighted_median",]$Var1 <- "KMA (TPM)"
dfp[dfp$Var2 == "kma_tpm_weighted_median",]$Var2 <- "KMA (TPM)"
dfp[dfp$Var1 == "kma_retention_weighted_median",]$Var1 <- "KMA (IR)"
dfp[dfp$Var2 == "kma_retention_weighted_median",]$Var2 <- "KMA (IR)"
dfp[dfp$Var1 == "irfinders_coverage_weighted_median",]$Var1 <- "IRFinder-S (coverage)"
dfp[dfp$Var2 == "irfinders_coverage_weighted_median",]$Var2 <- "IRFinder-S (coverage)"
dfp[dfp$Var1 == "irfinders_irratio_weighted_median",]$Var1 <- "IRFinder-S (IRratio)"
dfp[dfp$Var2 == "irfinders_irratio_weighted_median",]$Var2 <- "IRFinder-S (IRratio)"

lvlv <- c("iREAD (FPKM)", "IntEREst (IntRet)", "superintronic (score)", "KMA (TPM)",
          "KMA (IR)", "IRFinder-S (coverage)", "IRFinder-S (IRratio)")
dfp$Var1 <- factor(dfp$Var1, levels = lvlv)
dfp$Var2 <- factor(dfp$Var2, levels = lvlv)
# remove blank row and col
dfp <- dfp[!dfp$Var2 == "IRFinder-S (IRratio)",]
dfp <- dfp[!dfp$Var1 == "iREAD (FPKM)",]

# make the tile plot
ggtile <- ggplot(dfp, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black", na.rm = T) + theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.1,1), space = "Lab", 
                       name="Rho", na.value = "white") +
  theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
  geom_text(aes(label = label)) + xlab("") + ylab("") + ggtitle("iPSC")

ggtile

#-------------------------
# precision/recall results
#-------------------------
pr.fname <- "precision_recall_df_11-30-2021_12.46.22.csv"
pr <- data.table::fread(pr.fname, sep = ",", data.table = F)

# get threshold plot series
threshv <- c(0, 0.5, 0.9)
dfp <- do.call(rbind, lapply(threshv, function(ti){
  prf <- pr[pr$threshold==ti, grepl("precision|recall", colnames(pr))] # get p/r at thresh=ti
  dfp <- data.frame(precision = as.numeric(prf[grepl("precision", names(prf))]),
                    recall = as.numeric(prf[grepl("recall", names(prf))]),
                    algorithm = as.character(unique(gsub("_.*", "", names(prf)))),
                    stringsAsFactors = F);dfp$thresh <- ti; return(dfp)}))

ptplot <- ggplot(dfp, aes(x = dfp$precision, y = dfp$recall)) +
  geom_point(aes(color = algorithm)) + theme_bw() + 
  xlab("Precision") + ylab("Recall")

ptplot + facet_grid(cols = vars(thresh))

# get single plot at thresh = 0
threshi = 0
ptplot <- ggplot(dfp[dfp$thresh==threshi,], aes(x = precision, y = recall)) +
  geom_point(aes(color = algorithm)) + theme_bw() + 
  xlab("Precision") + ylab("Recall") + ggtitle(paste0("Threshold = ", threshi))
ptplot


#--------------
# distributions
#--------------
# z-score normalization
cnv <- colnames(irv@elementMetadata)
cnv.medians <- cnv[grepl("median$", cnv)]
dfp <- do.call(rbind, lapply(cnv.medians,function(cni){
  dat <- irv@elementMetadata[,cni]
  meani <- mean(dat, na.rm = T)
  sdi <- sd(dat, na.rm = T)
  dat.norm <- (dat-meani)/sdi
  dfi <- data.frame(dat = dat, dat.norm = dat.norm,
                    metric = rep(cni, length(dat)),
                    stringsAsFactors = F)
  return(dfi)}))
dfp$dat <- as.numeric(dfp$dat)
dfp$dat.norm <- as.numeric(dfp$dat.norm)

# violin plots -- whole dist
ggplot(dfp, aes(x = metric, y = dat.norm)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw()

# violin plots -- zoom
ggplot(dfp, aes(x = metric, y = dat.norm)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  ylim(0, 10)

#--------------
# recurrent top -- 99th percentile
#--------------
# 99th percentile filter
cnv <- colnames(irv@elementMetadata)
cnv.medians <- cnv[grepl("median$", cnv)]
lset <- lapply(cnv.medians, function(ci){
  dat <- irv@elementMetadata[,ci]
  qfilt <- quantile(dat, seq(0,1,0.01), na.rm = T)[100]
  which.ranges <- which(dat >= qfilt)
  return(irv[which.ranges])
}); names(lset) <- cnv.medians

length(subsetByOverlaps(lset[[1]], lset[[2]])) # 25
length(subsetByOverlaps(lset[[1]], lset[[3]])) # 1
length(subsetByOverlaps(lset[[1]], lset[[4]])) # 19
length(subsetByOverlaps(lset[[2]], lset[[3]])) # 26
length(subsetByOverlaps(lset[[2]], lset[[4]])) # 87
length(subsetByOverlaps(lset[[3]], lset[[4]])) # 21

# upset plot
lset.upset <- lapply(seq(length(lset)), function(li){
  paste0(seqnames(lset[[li]]), ":", 
         start(lset[[li]]), "-", end(lset[[li]]))
}); names(lset.upset) <- names(lset)
upset(fromList(lset.upset))

pdf("upset_99th-quant_ipsc.pdf", 5, 2)
upset(fromList(lset.upset)); dev.off()

# correlations among filtered introns
dati <- as.data.frame(irv, stringsAsFactors = F)
dati$intron.id <- paste0(dati$seqnames, ":", dati$start, "-", dati$end)
datf <- dati[dati$intron.id %in% unlist(lset.upset),]
dim(datf) # [1] 1382   20
datf <- datf[!is.na(datf$fract.expr.median),]
dim(datf) # [1] 1329   20
cor(datf[,grepl("median$", colnames(datf))],method = "spearman")

# recurrent ranges -- interest, superintronic, kma
recurr.str <- intersect(lset.upset[[2]], 
                        intersect(lset.upset[[3]], 
                                  lset.upset[[4]]))

"chr3:184182562-184183788" %in% recurr.str

#-----------
# bam pileup
#-----------

library(data.table); library(Gviz); library(rtracklayer)
library(Rsamtools); library(GenomicFeatures)
library(GenomicAlignments)

# load all transcript-read data for long read sample
ts.fname <- "RI_txs_to_read_ids_final.tsv"
ts <- fread(ts.fname, sep = "\t", header = T, data.table = F)

# get query info
query.str <- recurr.str[1]
chr <- gsub(":.*", "", query.str)
start <- gsub("(.*:|-.*)", "", query.str)
end <- gsub(".*-", "", query.str)

# get the corresponding transcript ids
coord.str <- paste0(".*",start, "-", end,".*")
tsf <- ts[grepl(coord.str, ts$tx_introns) & ts$chrom == chr,]
tid <- unique(tsf$transcript)
# get transcript coords from gtf
gtf.fpath <- file.path("gencode38", "gencode.v38.annotation.gtf")
gtf <- rtracklayer::import(gtf.fpath)
# gtf.df <- as.data.frame(gtf)
gtff <- gtf[gtf$transcript_id %in% tid]
length(gtff) # [1] 33

# get gviz tracks
gtff.df <- as.data.frame(gtff, stringsAsFactors = F)
genome = "hg38"
stack.genes = "full"
size.gtrack = 1
size.itrack = 1
size.genes = 3
tracks <- list(); message("Getting region info from gtff.df...")
chr <- as.character(unique(gtff.df$seqnames)) # get the chr, seqname
coord.values <- c(gtff.df$start, gtff$end) # get window coords
transcript.start <- min(coord.values)
transcript.end <- max(coord.values)

#message("Building tracks for region ", chr, ":", start, "-", end,"...")
message("Getting the Gviz tracks...")
tracks[["main.str"]] <- paste0(chr, ":", transcript.start, "-", transcript.end)
tracks[["gtrack"]] <- GenomeAxisTrack()
tracks[["itrack"]] <- IdeogramTrack(genome = genome, chromosome = chr)

# make the full transcript plot
fname = "newplot.pdf"
width = 8
height = 20
byrun = T
size.bamtrack = NULL
genome = "hg38"
tidv = paste0(tid, collapse = ";")
size.gtrack = 1
size.itrack = 1
size.genes = 3
sr.bam.fname <- "SRR6026510.sorted.bam"
lr.bam.fname <- "SRP098984_SAMN07611993.merged.aligned.sorted.bam"

plot.start <- as.numeric(transcript.start)
plot.end <- as.numeric(transcript.end)
tracks[["genes"]] <- GeneRegionTrack(gtff.df, genome = genome, 
                                     chromosome = chr, 
                                     name = "Gene Model", 
                                     start = plot.start, 
                                     end = plot.end,
                                     transcriptAnnotation = "transcript",
                                     stacking = "dense")
tracks[["atrack_bam_sr"]] <- AlignmentsTrack(sr.bam.fname, isPaired = TRUE, 
                                             genome = genome,
                                             start = plot.start, 
                                             end = plot.end, 
                                             chromosome = chr, 
                                             stacking = "dense")
tracks[["atrack_bam_lr"]] <- AlignmentsTrack(lr.bam.fname, isPaired = TRUE,
                                             genome = genome,
                                             start = plot.start,
                                             end = plot.end,
                                             chromosome = chr,
                                             stacking = "dense")
# get intron highlight
which.tracks <- which(names(tracks) %in% c("genes", "atrack_bam_sr", "atrack_bam_lr"))
ht <- HighlightTrack(trackList = tracks[which.tracks],
                     start = min(as.numeric(c(start, end))), 
                     width = abs(as.numeric(end)-as.numeric(start)),
                     chromosome = chr)
options(ucscChromosomeNames=FALSE)
pdf(fname, 5, 5)
plotTracks(ht)
dev.off()

# make the region plot
# make the full transcript plot
fname = "newplot_region.pdf"
width = 8
height = 20
byrun = T
size.bamtrack = NULL
genome = "hg38"
tidv = paste0(tid, collapse = ";")
size.gtrack = 1
size.itrack = 1
size.genes = 3
sr.bam.fname <- "SRR6026510.sorted.bam"
lr.bam.fname <- "SRP098984_SAMN07611993.merged.aligned.sorted.bam"

tracks <- list()
plot.start <- as.numeric(start)
plot.end <- as.numeric(end)
tracks[["genes"]] <- GeneRegionTrack(gtff.df, genome = genome, 
                                     chromosome = chr, 
                                     name = "Gene Model", 
                                     start = plot.start, 
                                     end = plot.end,
                                     transcriptAnnotation = "transcript",
                                     stacking = "dense")
tracks[["atrack_bam_sr"]] <- AlignmentsTrack(sr.bam.fname, isPaired = TRUE, 
                                             genome = genome,
                                             start = plot.start, 
                                             end = plot.end, 
                                             chromosome = chr, 
                                             stacking = "dense")
tracks[["atrack_bam_lr"]] <- AlignmentsTrack(lr.bam.fname, isPaired = TRUE,
                                             genome = genome,
                                             start = plot.start,
                                             end = plot.end,
                                             chromosome = chr,
                                             stacking = "dense")
# get intron highlight
options(ucscChromosomeNames=FALSE)
pdf(fname, 5, 5)
plotTracks(c(tracks[["genes"]], 
             tracks["atrack_bam_sr"], 
             tracks["atrack_bam_lr"]))
dev.off()





