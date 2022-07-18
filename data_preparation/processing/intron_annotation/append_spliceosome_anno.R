#!/usr/bin/env R

# Author: Sean Maden
#
# Append U2/U12 annotations to main harmonized tables.

library(GenomicRanges)

#----------
# load data
#----------
# get anno
u12 <- read.table("GRCh38_U12.bed", sep = "\t")
u2 <- read.table("GRCh38_U2.bed", sep = "\t")

# get intron tables
plot.titlev <- c("HX1", "iPSC")
tsv.fname.ipsc <- "called_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv"
tsv.fname.hx1 <- "called_RI_data_summary_HX1featureannotated_GCcontent.tsv"
ltsv <- list()
ltsv[["iPSC"]] <- read.table(tsv.fname.ipsc, sep = "\t", header = T)
ltsv[["HX1"]] <- read.table(tsv.fname.hx1, sep = "\t", header = T)

#------------------
# harmonize granges
#------------------
# get granges
# anno
colnames(u2) <- colnames(u12) <- c("chr", "start", "end", "feature", "v5", "strand")
u2$chr <- paste0("chr", u2$chr); u12$chr <- paste0("chr", u12$chr); 
u2gr <- makeGRangesFromDataFrame(u2, keep.extra.columns = T)
u12gr <- makeGRangesFromDataFrame(u12, keep.extra.columns = T)

# harmonize granges
tsv.hx1 <- ltsv$HX1
intronv <- tsv.hx1$intron
hx1gr <- makeGRangesFromDataFrame(data.frame(chr = gsub(":.*", "", intronv), 
                                             start = gsub(".*:|-.*", "", intronv),
                                             end = gsub(".*-", "", intronv)))
tsv.ipsc <- ltsv$iPSC
intronv <- tsv.ipsc$intron
ipscgr <- makeGRangesFromDataFrame(data.frame(chr = gsub(":.*", "", intronv), 
                                             start = gsub(".*:|-.*", "", intronv),
                                             end = gsub(".*-", "", intronv)))

#------------
# append anno
#------------
# hx1
u2.fol <- findOverlaps(hx1gr, u2gr)
u12.fol <- findOverlaps(hx1gr, u12gr)
annov <- is.u2 <- is.u12 <- rep("other", length(hx1gr))
annov[queryHits(u2.fol)] <- "u2"; annov[queryHits(u12.fol)] <- "u12"
annov[intersect(queryHits(u2.fol),queryHits(u12.fol))] <- "other"
tsv.hx1$intron_type_annotation <- annov
# ipsc
u2.fol <- findOverlaps(ipscgr, u2gr)
u12.fol <- findOverlaps(ipscgr, u12gr)
annov <- is.u2 <- is.u12 <- rep("other", length(ipscgr))
annov[queryHits(u2.fol)] <- "u2"; annov[queryHits(u12.fol)] <- "u12"
annov[intersect(queryHits(u2.fol),queryHits(u12.fol))] <- "other"
tsv.ipsc$intron_type_annotation <- annov

# save new files
hx1.cont.fname = "called_RI_data_summary_HX1featureannotated_GCcontent_splicetype.tsv"
write.table(tsv.hx1, file = hx1.cont.fname, row.names = F, sep = "\t")
# ipsc
ipsc.cont.fname = "called_RI_data_summary_iPSCfeatureannotated_GCcontent_splicetype.tsv"
write.table(tsv.ipsc, file = ipsc.cont.fname, row.names = F, sep = "\t")
