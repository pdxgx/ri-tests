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
srrid <- "SRR2911306"
run.handle <- "hx1"

# load the majiq results
mfname <- paste0(run.handle, "_majiq.psi.tsv")
mfdir <- paste0(srrid, "_", run.handle)
mfpath <- file.path(mfdir, mfname)
mj <- read.table(mfpath, sep = "\t", header = T)

# subset ri
mj <- mj[!mj$ir_coords=="",]
dim(mj) # [1] 12552     9

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
dim(mj) # [1] 12552    11
mj <- mj[!mj$chr=="NA",]
dim(mj) # [1] 12537    11
# summary
table(mj$chr)
# 1         10         11         12         13         14         15         16         17         18         19          2 
# 1386        361        843        711        182        412        394        800       1002         88        997        675 
# 20         21         22          3          4          5          6          7          8          9 GL000219.1          X 
# 344        155        408        632        304        530        637        573        371        400          1        321 
# Y 
# 10
mj[c(1:10),c(1,10,11)]
# gene_id chr        gene.new
# 2  ENSG00000125637.16   2 ENSG00000125637
# 4  ENSG00000125637.16   2 ENSG00000125637
# 8  ENSG00000198836.10   3 ENSG00000198836
# 9  ENSG00000112531.17   6 ENSG00000112531
# 11 ENSG00000112531.17   6 ENSG00000112531
# 16 ENSG00000112531.17   6 ENSG00000112531
# 25 ENSG00000147535.17   8 ENSG00000147535
# 26  ENSG00000254402.7   8 ENSG00000254402
# 29 ENSG00000125637.16   2 ENSG00000125637
# 34  ENSG00000254402.7   8 ENSG00000254402

#---------------------------
# collapse repeated features
#---------------------------
# note: use median psi across repeated features
# get ri stats
mj$majiq.ir.mean.psi <- gsub(".*;", "", mj$mean_psi_per_lsv_junction)
# mj$ir.sd.psi <- gsub(".*;", "", mj$stdev_psi_per_lsv_junction)e

# check duplicates
dft <- as.data.frame(table(mj$ir_coords))
length(unique(dft[,1])) # 11876
ir.coord.dupv <- unique(dft[dft[,2]>1,1])
length(ir.coord.dupv) # [1] 647

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
dim(dfir) # [1] 11876     3
head(dfir)
# ir_coord chr majiq.ir.mean.psi
# 1 113193949-113195726   2            0.0250
# 2 113196608-113197194   2            0.0294
# 3 193646580-193647064   3            0.2009
# 4 163563720-163564666   6            0.1353
# 5 163566796-163569437   6            0.3903
# 6 163564711-163565945   6            0.1644

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