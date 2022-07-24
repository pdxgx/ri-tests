#!/usr/bin/env R

# Author: Sean Maden
#
# Test for circRNA overlap with TP, FP, FN intron truth categories.

library(GenomicRanges)

#----------
# load data
#----------
# load intron truth category data
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "called_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv"
tsv.fname.hx1 <- "nonzero_RI_data_summary_HX1featureannotated_GCcontent"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

# get cRNA ranges
crna <- read.table("intronic_circRNAs_GRCh38.txt", sep = " ")
crna$chr <- gsub(":.*", "", crna[,1])
crna$start <- gsub(".*:|-.*", "", crna[,1])
crna$end <- gsub(".*-", "", crna[,1])
crna.gr <- makeGRangesFromDataFrame(crna)

# color palettes
pal <- c("FN" = "#6d9cc6", "FP" = "#d8788a", "TP" = "#9db92c")
pal.sr <- c('IRFinder-S' = '#e1665d', 'superintronic' = '#f8b712', 
         'iREAD' = '#689404', 'IntEREst' = '#745bad', 'KMA' = '#33a2b7',
         'rMATS' = '#DAF7A6', 'MAJIQ' = '#FFC0C0', 'SUPPA2' = '#abddff')

#-----------------
# helper functions
#-----------------
# get intronv by group label
get_singlegroup <- function(tsv, min.tmetric = 4, regex = "true_positives$"){
  tsvf <- tsv[!duplicated(tsv$intron),]
  tsvf <- tsv[,c(colnames(tsv)[grepl(regex, colnames(tsv))], "intron")]
  group.data <- tsvf[apply(tsvf, 1, function(ri){
    length(ri[ri==1])>=min.tmetric}),]$intron
  return(group.data)
}

# get intronv groups 
get_lgroup <- function(tsv, min.tmetric = 4,
                       regexv = c("true_positives$", "false_negatives$", "false_positives$")){
  lgroup <- list()
  grp.namev <- paste0(gsub("\\$", "", regexv), "_", min.tmetric)
  lgroup <- lapply(regexv, function(regex){
    get_singlegroup(tsv, regex, min.tmetric = min.tmetric)
  })
  names(lgroup) <- grp.namev
  return(lgroup)
}

#------------------------------
# get crna overlap -- by srtool
#------------------------------
toolv <- c("iread", "interest", "superintronic", "kma", "irfinders", 
           "majiq", "rmats", "suppa2")

# get sets
# get filtered intron ids/coords
lgroup <- lapply(ltsv, function(tsvi){
  lgroupi <- lapply(toolv, function(tooli){
    tsvf <- tsvi[!duplicated(tsvi$intron),]
    filt.col <- grepl(".*filtintron_lwm$", colnames(tsvf)) & 
      grepl(tooli, colnames(tsvf))
    which.nz <- which(tsvf[,filt.col] > 0 & !is.na(tsvf[,filt.col]))
    return(tsvf[which.nz,]$intron)
  })
  names(lgroupi) <- toolv
  lgroupi
})
names(lgroup) <- names(ltsv)

# parse overlaps by intron coords
dfol.sr <- do.call(rbind, lapply(seq(length(lgroup)), function(ii){
  samplei <- names(lgroup)[ii]; lgroupi <- lgroup[[samplei]]
  do.call(rbind, lapply(names(lgroupi), function(namei){
    # get the grset
    intronv <- lgroupi[[namei]]
    dfintv <- data.frame(start = gsub(".*:|-.*", "", intronv),
                         end = gsub(".*-", "", intronv),
                         chr = gsub(":.*", "", intronv))
    intronv.gr <- makeGRangesFromDataFrame(dfintv)
    # get df of overlap summaries
    num.intron.ol <- length(subsetByOverlaps(intronv.gr, crna.gr))
    num.crna.ol <- length(subsetByOverlaps(crna.gr, intronv.gr))
    data.frame(tmetric = namei, sample = samplei, total.intron = length(intronv.gr),
               num.intron.ol = num.intron.ol, fract.intron.ol = num.intron.ol/length(intronv.gr),
               num.crna.ol = num.crna.ol)
  }))
}))

# make new barplots
# get plot data
dfpr <- dfol.sr
lvlv <- c("IntEREst", "KMA", "iREAD", "superintronic", "IRFinder-S",
          "MAJIQ", "rMATS", "SUPPA2")
dfpr$Tool <- factor(dfpr$tmetric, levels = lvlv)
dfpr[dfpr$tmetric=="interest",]$Tool <- "IntEREst"
dfpr[dfpr$tmetric=="iread",]$Tool <- "iREAD"
dfpr[dfpr$tmetric=="kma",]$Tool <- "KMA"
dfpr[dfpr$tmetric=="irfinders",]$Tool <- "IRFinder-S"
dfpr[dfpr$tmetric=="majiq",]$Tool <- "MAJIQ"
dfpr[dfpr$tmetric=="rmats",]$Tool <- "rMATS"
dfpr[dfpr$tmetric=="suppa2",]$Tool <- "SUPPA2"
# format plot data
dfpr$perc.intron.ol <- 100*dfpr$fract.intron.ol
dfpr$`Truth\nmetric` <- dfpr$tmetric.label
# get new plot objects
gg.crna <- ggplot(dfpr, aes(x = Tool, y = perc.intron.ol, fill = Tool)) +
  geom_bar(stat = "identity") + theme_bw() + 
  ylab("Percent of introns\noverlapping cRNAs") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = pal.sr)
gg.crna <- gg.crna + facet_wrap(~sample, nrow = 2)

# save new figures
plot.fname <- "ggbp_crna-overlaps_by-srtool_combined-2samples"
width <- 2; height <- 2.5
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, width, height)
print(gg.crna); dev.off()
# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = width, height = height, res = 500, units = "in")
gg.crna; dev.off()

#-------------------------------------------------
# get crna overlap -- by intron truth metric group
#-------------------------------------------------
# get sets
dfol <- do.call(rbind, lapply(seq(5), function(si){
  lgroup <- lapply(ltsv, function(tsvi){get_lgroup(tsvi, si)})
  names(lgroup) <- names(ltsv)
  do.call(rbind, lapply(seq(length(lgroup)), function(ii){
    samplei <- names(lgroup)[ii]; lgroupi <- lgroup[[samplei]]
    do.call(rbind, lapply(names(lgroupi), function(namei){
      # get the grset
      intronv <- lgroupi[[namei]]
      dfintv <- data.frame(start = gsub(".*:|-.*", "", intronv),
                           end = gsub(".*-", "", intronv),
                           chr = gsub(":.*", "", intronv))
      intronv.gr <- makeGRangesFromDataFrame(dfintv)
      # get df of overlap summaries
      num.intron.ol <- length(subsetByOverlaps(intronv.gr, crna.gr))
      num.crna.ol <- length(subsetByOverlaps(crna.gr, intronv.gr))
      data.frame(tmetric = namei, sample = samplei, total.intron = length(intronv.gr),
                 num.intron.ol = num.intron.ol, fract.intron.ol = num.intron.ol/length(intronv.gr),
                 num.crna.ol = num.crna.ol)
    }))
  }))
}))

# save new data, table
dfol.fname <- "dfol-crna_tmetric-varfilt_combined-2samples"
write.csv(dfol, file = paste0(dfol.fname, ".csv"))
save(dfol, file = paste0(dfol.fname, "rda"))

#------------------
# save new boxplots
#------------------
# color palette
pal <- c("FN" = "#6d9cc6", "FP" = "#d8788a", "TP" = "#9db92c")

# format plot data
dfol$perc.intron.ol <- 100*dfol$fract.intron.ol
dfol$tmetric.label <- ifelse(grepl("true_positives", dfol$tmetric), "TP",
                             ifelse(grepl("false_negatives", dfol$tmetric), "FN", "FP"))
dfol$tmetric.label <- factor(dfol$tmetric.label, levels = c("FP", "FN", "TP"))
dfol$`Truth\nmetric` <- dfol$tmetric.label
# get new plot objects
gg.crna <- ggplot(dfol, aes(x = tmetric.label, y = perc.intron.ol, 
                            fill = `Truth\nmetric`)) +
  geom_boxplot() + theme_bw() + ylab("Percent of introns\noverlapping cRNAs") +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = pal)
gg.crna <- gg.crna + facet_wrap(~sample, nrow = 2)

# save new figures
plot.fname <- "ggbp_crna-fract-ol_combined-2samples"
# new pdf
pdf.fname <- paste0(plot.fname, ".pdf")
pdf(pdf.fname, 4, 3)
print(gg.crna); dev.off()
# new png
png.fname <- paste0(plot.fname, ".png")
png(png.fname, width = 4, height = 3, res = 500, units = "in")
gg.crna; dev.off()