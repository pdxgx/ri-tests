#!/usr/bin/env R

# Author: Sean Maden
#
# Plot the reads mapped by run/sample (STAR and bowtie2)

library(ggplot2)
library(ggpubr)
library(gridExtra)

#-----------------
# helper functions
#-----------------

# make composite barplot
bpcomp <- function(dfp, legend.title = "Value", ylab1.str = "Count", ylab2.str = "Percent",
                   pdf.fname = "new.pdf", pdf.width = 5, pdf.height = 3,
                   ylab.str = paste0("Sample", paste0(rep(" ", 25), collapse = ""), 
                                      collapse = "")){
  # get plot objects
  bp.count <- ggplot(dfp, aes(x = sample, y = number, fill = value)) +
    geom_bar(stat = "identity", colour = "black") + theme_bw() + ylab(ylab1.str) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  bp.perc <- ggplot(dfp, aes(x = sample, y = number, fill = value)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") + 
    theme_bw() + ylab(ylab2.str) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  plot.legend <- ggplot(dfp, aes(x = sample, y = number, fill = value)) +
    geom_bar(stat = "identity") + guides(fill=guide_legend(title=legend.title))
  bp.legend <- get_legend(plot.legend)
  
  # save new plot
  pdf(pdf.fname, width = pdf.width, height = pdf.height)
  grid.arrange(bp.count, bp.perc, bp.legend, nrow = 1, bottom = ylab.str)
  dev.off()
  return(list(bp.count = bp.count, bp.perc = bp.perc, 
              bp.legend = bp.legend))
}

#------------------------------
# long read lima report summary
#------------------------------
# read in line data
dpath <- "lima_reports"; lr <- list()

dnv <- list.files(dpath)
lr <- lapply(dnv, function(type){
  fnv <- list.files(file.path(dpath, type))
  li <- lapply(fnv, function(file){
    rfpath <- file.path(dpath, type, file)
    readLines(rfpath)}); names(li) <- fnv; return(li)}); names(lr) <- dnv

for(type in list.files(dpath)){
  lr[[type]] <- list()
  for(file in list.files(file.path(dpath, type))){
    rfpath <- file.path(dpath, type, file)
    lr[[type]][[file]] <- readLines(rfpath)
  }
}
# format data
dfp <- do.call(rbind, lapply(names(lr), function(typei){lri <- lr[[typei]]
  do.call(rbind, lapply(names(lri), function(reporti){
    lrii <- lri[[reporti]]
    dfi <- data.frame(total = as.numeric(gsub(".*:| ", "", lrii[1])),
               above.thresh = as.numeric(gsub(".*:|\\(.*| ", "", lrii[2])),
               below.thresh.minlen = as.numeric(gsub(".*:|\\(.*| ", "", lrii[6])),
               below.thresh.minscore = as.numeric(gsub(".*:|\\(.*| ", "", lrii[7])),
               below.thresh.minendscore = as.numeric(gsub(".*:|\\(.*| ", "", lrii[8])),
               below.thresh.minpasses = as.numeric(gsub(".*:|\\(.*| ", "", lrii[9])),
               below.thresh.minscorelead = as.numeric(gsub(".*:|\\(.*| ", "", lrii[10])),
               below.thresh.minrefspan = as.numeric(gsub(".*:|\\(.*| ", "", lrii[11])),
               below.thresh.withoutsmrtbelladapter = as.numeric(gsub(".*:|\\(.*| ", "", lrii[12])),
               below.thresh.undesired5p5ppairs = as.numeric(gsub(".*:|\\(.*| ", "", lrii[13])),
               below.thresh.undesired3p3ppairs = as.numeric(gsub(".*:|\\(.*| ", "", lrii[14])),
               stringsAsFactors = F)
    dfi$run.handle <- typei; dfi$srrid <- gsub("_.*|\\..*", "", reporti)
    return(dfi)}))}))

# barplots, threshold medians
thresh.tot.med.ipsc <- median(dfp[dfp$run.handle=="ipsc",]$total)
thresh.pos.med.ipsc <- median(dfp[dfp$run.handle=="ipsc",]$above.thresh)
thresh.tot.med.hx1 <- median(dfp[dfp$run.handle=="hx1",]$total)
thresh.pos.med.hx1 <- median(dfp[dfp$run.handle=="hx1",]$above.thresh)
numv <- c(thresh.pos.med.ipsc, thresh.tot.med.ipsc - thresh.pos.med.ipsc,
          thresh.pos.med.hx1, thresh.tot.med.hx1 - thresh.pos.med.hx1)
dfpm <- data.frame(sample = rep(c("ipsc", "hx1"), each = 2), value = rep(c(1,0), 2), 
           number = numv, stringsAsFactors = F)
dfpm$value <- factor(dfpm$value)
# plot threshold medians
bpcomp(dfpm, legend.title = "Above lima\nthreshold", 
       ylab1.str = "Read count", ylab2.str = "Read percent", 
       pdf.fname = "A_bp-lima-reads-thresh_2samples.pdf")

# barplots, median flag counts, below thresholds
dff <- do.call(rbind, lapply(unique(dfp$run.handle), function(ri){
  dfp.filt <- dfp$run.handle == ri; mdat <- dfp[dfp.filt,c(3:11)]
  if(ri == "ipsc"){mdat <- t(apply(mdat,2,median))}
  mdat <- as.data.frame(t(mdat), stringsAsFactors = F)
  mdat$flag <- rownames(mdat); mdat$run.handle <- ri; 
  colnames(mdat) <- c("num.reads", "flag", "run.handle")
  return(mdat)}))
# format plot data
dff <- dff[!dff[,1]==0,]
dff$flag <- gsub("below.thresh.", "", dff$flag)
lvlv <- c("minendscore", "undesired3p3ppairs", "minrefspan", "undesired5p5ppairs", "minlen")
dff$flag <- factor(dff$flag, levels = lvlv)
# get plot object
bp <- ggplot(dff, aes(x = run.handle, y = num.reads, fill = flag)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
  guides(fill=guide_legend(title="Below lima\nthresh. flag")) +
  ylab("Read count") + xlab("Sample")
# save plot
pdf.fname <- "B_bp-lima-subthresh-flags_2samples.pdf"
pdf(pdf.fname, 5, 3); print(bp); dev.off()

#---------------------
# star -- reads mapped
#---------------------
# get plot data
dfp <- data.frame(sample = c("ipsc", "ipsc", "hx1", "hx1"),
           value = factor(rep(c(0, 1), 2)),
           number = c(91330785 - 53916311, 53916311,
                         24463210 - 21618226, 21618226),
           stringsAsFactors =F)

# get plot
lbp <- bpcomp(dfp, legend.title = "Uniquely\nmapping", 
              ylab1.str = "Read count", ylab2.str = "Read percent",
              pdf.fname = "C_bp-star-reads-umap_2samples.pdf")

#----------------------------
# star -- junctions annotated
#----------------------------
# get plot data
dfjx <- data.frame(sample = c(rep(c("ipsc", "hx1"), each = 2)),
           value = factor(c(rep(c(0,1), 2))),
           number = c(105593, 202939, 52885, 162992), 
           stringsAsFactors = F)

# get plot
lbp <- bpcomp(dfjx, legend.title = "Annotated", 
              ylab1.str = "Jx count", ylab2.str = "Jx percent",
              pdf.fname = "D_bp-star-jx-anno_2samples.pdf")

#------------------------
# bowtie2 -- reads mapped
#------------------------
# note: reads mapped to intron regions for KMA pipelilne
dfbt2 <- data.frame(sample = c(rep("ipsc", 3), rep("hx1", 3)),
                      value = factor(rep(c("0", "1", ">1"),2), 
                                     levels = c("0", "1", ">1")),
                      number = c(41470026, 12388302, 37472457,
                                    3441596, 4934735, 16086879),
                      stringsAsFactors = F)

# get plot
lbp <- bpcomp(dfbt2, legend.title = "Times\naligned", 
              ylab1.str = "Read count", ylab2.str = "Read percent",
              pdf.fname = "E_bp-bowtie2-read-timesaligned_2samples.pdf")

#----------------
# irfinders flags
#----------------
# load data
hx1.fpath <- file.path("SRR2911306_hx1",
                       "df-granges_irfinders_SRR2911306-hx1.csv")
ipsc.fpath <- file.path("SRR6026510_ipsc",
                        "df-granges_irfinders_SRR6026510-ipsc.csv")
hx1 <- data.table::fread(hx1.fpath, data.table = F, sep = ",")
ipsc <- data.table::fread(ipsc.fpath, data.table = F, sep = ",")
lsamp <- list(hx1 = hx1, ipsc = ipsc)

# get barplot dfp
warningsv <- unique(hx1$Warnings)
dfflag <- do.call(rbind, lapply(warningsv, function(flag){
  do.call(rbind, lapply(names(lsamp), function(sample){
    dat <- lsamp[[sample]]
    data.frame(value = flag, number = nrow(dat[dat$Warnings==flag,]), 
               sample = sample, stringsAsFactors = F)}))}))
# format vars
dfflag[dfflag$value == "-",]$value <- "None"
lvlv <- c("LowSplicing", "MinorIsoform", "NonUniformIntronCover", 
          "None", "LowCover")
dfflag$value <- factor(dfflag$value, levels = lvlv)

# get plot object
bp <- ggplot(dfflag, aes(x = sample, y = number, fill = value)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() + 
  xlab("Sample") + ylab("Intron count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
  guides(fill=guide_legend(title="IRFinder-S\nflag"))
# save plot
pdf.fname <- "F_bp-irfinders-intronflags_2samples.pdf"
pdf(pdf.fname, 5, 3); print(bp); dev.off()
