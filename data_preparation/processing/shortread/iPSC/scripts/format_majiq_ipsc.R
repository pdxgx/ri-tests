#!/usr/bin/env R

# Author: Sean Maden
# 
# Format the MAJIQ results and store as GenomicRanges object.
#

library(GenomicRanges)
library(biomaRt)

#----------
# load data
#----------
# set the run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

# load the majiq results
mfname <- paste0(run.handle, "_majiq.psi.tsv")
mfdir <- paste0(srrid, "_", run.handle)
mfpath <- file.path(mfdir, mfname)
mj <- read.table(mfpath, sep = "\t", header = T)

# subset ri
mj <- mj[!mj$ir_coords=="",]
dim(mj) # [1] 15979     9

# load mart file
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#------------------------
# map data to chromosomes
#------------------------
# get chrs for unique genes
genev <- unique(gsub("\\..*", "", mj$gene_id))
length(genev) # 5399
chrv <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name'), 
              filters = c('ensembl_gene_id'), values = genev, mart = ensembl)

# append chr
mj$gene.new <- gsub("\\..*", "", mj$gene_id)
mj$chr <- sapply(seq(nrow(mj)), function(ii){
  ifelse(mj$gene.new[ii] %in% chrv[,1],
         chrv[chrv[,1]==mj$gene.new[ii],2], "NA")
})

# filter missing chr
dim(mj) # [1] 15979    11
mj <- mj[!mj$chr=="NA",]
dim(mj) # [1] 15957    11

# summaries
table(mj$chr)
# 1         10         11         12         13         14         15         16         17         18         19          2 
# 1561        466       1016        861        220        522        577        980       1139        201       1405        857 
# 20         21         22          3          4          5          6          7          8          9 GL000195.1 KI270711.1 
 #598        229        520        774        438        685        638        747        518        594          1          2 
# KI270721.1 KI270734.1          X 
# 1          7        400
mj[c(1:10),c(1,10,11)]
#               gene_id        gene.new chr
# 12 ENSG00000198064.13 ENSG00000198064  16
# 19 ENSG00000077454.16 ENSG00000077454   7
# 26 ENSG00000198064.13 ENSG00000198064  16
# 30 ENSG00000184402.15 ENSG00000184402  20
# 38 ENSG00000204348.10 ENSG00000204348   6
# 41 ENSG00000184402.15 ENSG00000184402  20
# 43 ENSG00000183506.17 ENSG00000183506  22
# 59 ENSG00000145740.19 ENSG00000145740   5
# 61 ENSG00000184402.15 ENSG00000184402  20
# 62 ENSG00000184402.15 ENSG00000184402  20

#---------------------------
# collapse repeated features
#---------------------------
# note: use median psi across repeated features
# get ri stats
mj$majiq.ir.mean.psi <- gsub(".*;", "", mj$mean_psi_per_lsv_junction)
# mj$ir.sd.psi <- gsub(".*;", "", mj$stdev_psi_per_lsv_junction)e

# check duplicates
dft <- as.data.frame(table(mj$ir_coords))
length(unique(dft[,1])) # [1] 14954
ir.coord.dupv <- unique(dft[dft[,2]>1,1])
length(ir.coord.dupv) # [1] 986

# check duplicate values
# mj[mj$ir_coords==ir.coord.dupv[1],]$ir.mean.psi

# make new dataset
dfir <- data.frame(ir_coord = mj$ir_coords, chr = mj$chr,
                   majiq.ir.mean.psi = mj$majiq.ir.mean.psi)
dfir <- dfir[!dfir$ir_coord %in% ir.coord.dupv,] # filter dfir
dim(dfir)
dfir <- rbind(dfir, do.call(rbind, lapply(ir.coord.dupv, function(ir.coordi){
  dfii <- mj[mj$ir_coords==ir.coordi,]
  medi <- median(as.numeric(dfii$majiq.ir.mean.psi), na.rm = T)
  data.frame(ir_coord = ir.coordi, chr = unique(dfii$chr),
             majiq.ir.mean.psi = medi)
})))
dim(dfir) # [1] 14954     3
head(dfir)
#               ir_coord chr majiq.ir.mean.psi
# 1   30245878-30254309  16            0.9583
# 2 100574357-100574659   7            0.9877
# 3   30240243-30245453  16            0.9300
# 5   31970379-31970441   6            0.9683
# 7   21475301-21475991  22            0.3815
# 8   69121753-69121753   5            0.9972

#--------------------
# make and save grset
#--------------------
# separate coords
dfir$start <- gsub("-.*", "", dfir$ir_coord)
dfir$stop <- gsub(".*-", "", dfir$ir_coord)
mj.gr <- makeGRangesFromDataFrame(dfir[,c(2:5)], keep.extra.columns = T)

# save granges
mj.gr.fname <- paste0("granges_majiq_", srrid,"-",run.handle,".rda")
mj.gr.fpath <- file.path(mfdir, mj.gr.fname)
save(mj.gr, file = mj.gr.fpath)
