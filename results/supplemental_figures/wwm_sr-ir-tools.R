#!/usr/bin/env R

# Author: Sean Maden
#
# Width-weighted medians figure

library(gridExtra)
library(ggplot2)

#----------
# load data
#----------
# load individual ranges
dpath <- "SRR2911306_hx1"
fnv <- c("granges_superintronic_SRR2911306-hx1.rda",
         "granges_irfinders_SRR2911306-hx1.rda",
         "granges_kma_SRR2911306-hx1.rda",
         "granges_iread_SRR2911306-hx1.rda",
         "granges_interest_SRR2911306-hx1.rda")
lgr <- lapply(fnv, function(fi){
  get(load(file.path(dpath, fi)))
})
names(lgr) <- gsub("^granges_|_SRR.*", "", fnv)

# get merged data
gr.fname <- "granges-lrmap_sr-5-methods_SRR2911306-hx1.rda"
grmerge <- get(load(file.path(dpath, gr.fname)))
em <- elementMetadata(grmerge)
em <- em[,grepl(".*lwm$", colnames(em))]
nz.filt <- apply(em, 1, function(ri){length(ri[ri>0])==7})
grmf <- grmerge[nz.filt]

#----------------------
# get plot data by tool
#----------------------
gr.sol <- findOverlaps(grmf, lgr[[4]])
qh.sol <- queryHits(gr.sol)
sh.sol <- subjectHits(gr.sol)

# get region with many overlaps
dft <- as.data.frame(table(qh.sol))
qi <- as.numeric(dft[which(dft[,2] > 2),1])

# get long read region/view window xaxis coords
lrgr.filt <- grmf[qi[8]]
min.range <- start(lrgr.filt)
max.range <- end(lrgr.filt)

# get vectors for tool iter 
tool.namev <- c("superintronic", "irfinders_irratio", "kma_tpm", "iread", "interest")
tool.labv <- c("superintronic", "IRFinder-S", 
               "KMA", "iREAD", "IntEREst")
colv <- c("purple", "red", "green", "blue", "orange")
metricv <- c("score", "IRratio", "kma.tpm", "fpkm", "IntRet_frequency")
metric.labelv <- c("Score", "IRratio", "TPM", "FPKM", "Freq.")

# get the wwm vector
em <- elementMetadata(lrgr.filt)
em <- em[,grepl(".*allintron.*", colnames(em))]
em <- em[,grepl(paste0("^",tool.namev,".*",collapse ="|"), colnames(em))]; em <- unlist(em)
em <- c(em[grepl(tool.namev[1], names(em))], em[grepl(tool.namev[2], names(em))],
        em[grepl(tool.namev[3], names(em))], em[grepl(tool.namev[4], names(em))],
        em[grepl(tool.namev[5], names(em))])
wwmv <- as.numeric(unlist(em)) # get WWM

#-----------------------
# make new plot objects
#-----------------------
lgg <- lapply(seq(length(lgr)), function(ii){
  message(ii)
  # tool iterables
  gri <- lgr[[ii]]; metrici <- metricv[ii]; wwm <- wwmv[ii]
  metric.labeli <- metric.labelv[ii]; coli <- colv[ii]
  tool.name <- tool.namev[ii]; tool.lab <- tool.labv[ii]
  # ggrect plot title str
  title.str <- paste0(tool.lab, " (LWM = ", round(wwm, 3), ")")
  gr.sol <- findOverlaps(lrgr.filt, gri)
  qh.sol <- queryHits(gr.sol); sh.sol <- subjectHits(gr.sol)
  # get the rectangle coords for the metric
  srgr.filt <- gri[sh.sol]
  dfp <- do.call(rbind, lapply(seq(length(srgr.filt)), function(jj){
    grj <- srgr.filt[jj]; coordv <- c(start(grj), end(grj))
    start.coord <- min(coordv); end.coord <- max(coordv)
    xmin <- ifelse(start.coord < min.range, min.range, start.coord)
    xmax <- ifelse(end.coord > max.range, max.range, end.coord)
    xmax <- ifelse(width(grj) == 1, xmin + 1, xmax)
    emj <- elementMetadata(grj)
    data.frame(xmin = xmin, xmax = xmax, ymin = 0, ymax = emj[,metrici])
  }))
  max.ylab <- max(dfp$ymax)
  ggrect <- ggplot() + geom_rect(data = dfp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                                 fill = "black", alpha = 0.5) +
    geom_hline(yintercept = as.numeric(wwm), alpha = 0.5, size = 0.5, color = "red") +
    xlim(min.range, max.range) + ggtitle(title.str) + theme_bw() +
    theme(axis.title.x = element_blank()) + ylab(metric.labeli) + ylim(0,max.ylab)
  
  pdf("newplot.pdf", 5, 3)
  ggsmooth <- ggplot(dfp, aes(x = xmin, y = ymax)) + geom_smooth()
  ggsmooth
  dev.off()
  
  return(list("rect" = ggrect))
})

#---------------
# save new plots
#---------------
# get plot data
gcoord.str <- paste0("Intron ",
                     as.character(seqnames(lrgr.filt)[1]), ":",
                     as.character(start(lrgr.filt)), "-",
                     as.character(end(lrgr.filt)))s
plot.fname <- "sfig_LWM-rect-coord-exe_sr-ir-5tools"
# save new pdf
pdf(paste0(plot.fname, ".pdf"), 4.5, 5.5)
grid.arrange(lgg[[1]][["rect"]], lgg[[2]][["rect"]], 
             lgg[[3]][["rect"]], lgg[[4]][["rect"]], 
             lgg[[5]][["rect"]], nrow = 5,
             bottom = "Coordinate (bp)", 
             top = gcoord.str)
dev.off()
# save new png
png(paste0(plot.fname, ".png"), width = 4.5, height = 5.5, res = 500, units = "in")
grid.arrange(lgg[[1]][["rect"]], lgg[[2]][["rect"]], 
             lgg[[3]][["rect"]], lgg[[4]][["rect"]], 
             lgg[[5]][["rect"]], nrow = 5,
             bottom = "Coordinate (bp)", 
             top = gcoord.str)
dev.off()