#!/usr/bin/env R

# Author: Sean Maden
# 
# Analyze the short read data, mapped to long read intron ranges.

library(ggplot2)
library(UpSetR)
library(tidyverse)

srrid <- "SRR2911306"; run.handle <- "hx1"

#----------
# load data
#----------
irv.fname <- paste0("granges-lrmap_sr-5-methods_", srrid,"-",run.handle,".rda")
irv <- get(load(irv.fname))

# get median gene coverage
mpile <- read.csv("mpile-sstat_gencode-v35_SRR2911306-hx1.csv")
grpile <- makeGRangesFromDataFrame(mpile, keep.extra.columns = T)
grpile <- grpile[!is.na(grpile$median)]
grpile <- grpile[grpile$median >= 2]
# filt irv on coverage
irvf <- subsetByOverlaps(irv, grpile)

#--------------
# get plot data
#--------------
# get metadata df
dati <- as.data.frame(irvf@elementMetadata, stringsAsFactors = F)

# get correlation matrix
lcor <- lapply(c("filtintron_wwm$", "allintron_wwm$"), function(filt){
  cor(dati[,which(grepl(filt, colnames(dati)))], method = "spearman")
}); names(lcor) <- c("filt", "all")

ggplot(lcor$filt, aes(x = ))

# get dfp list
ldfp <- lapply(seq(length(lcor)), function(ii){
  dfp <- lcor[[ii]]
  # append revised rownames
  rownames(dfp) <- c("iREAD (FPKM)", "IntEREst (FPKM)", "superintronic (score)",
                     "KMA (TPM)", "IRFinder-S (IRratio)")
  which.na <- which(upper.tri(dfp)|dfp==1); dfp[which.na] <- "NA"
  dsh <- reshape2::melt(dfp); dsh
}); names(ldfp) <- c("Filtered", "All")

lgg <- lapply(seq(length(ldfp)), function(ii){
  dfp <- ldfp[[ii]]
  dfp$value <- as.numeric(dfp$value)
  dfp$label <- as.character(round(dfp$value, 2))
  filt.na <- grepl("IRFinder-S", dfp$Var2)|grepl("iREAD", dfp$Var1)
  dfp <- dfp[!filt.na,]
  # make new plot object
  ggplot(dfp, aes(x = Var2, y = Var1, fill = value)) + geom_tile(color = "black") +
    theme_bw() + geom_text(aes(x = Var2, y = Var1, label = label), color = "white") +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle(names(ldfp)[ii])
})

#-----------------------
# save new corr heatmaps
#-----------------------


#-------------------------
# precision/recall results
#-------------------------
pr.fname <- "precision_recall_df_11-30-2021_12.47.27.csv"
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





