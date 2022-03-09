#!/usr/bin/env R

# Author: Sean Maden
#
# Run superintronic on BAM chr chunks of ips short read run.
#
# Remote server instructions:
# 1. run the conda env
# > conda activate 405_superinterest
#
# 2. run R with sudo access
# > sudo /opt/anaconda3/envs/405_superinterest/bin/R
#

library(superintronic)
library(plyranges)
library(Rsamtools)
#library(ggplot2)

# set run id
srrid <- "SRR2911306"

#----------
# load data
#----------
# manage paths
bamdpath <- file.path("eternity", "data", "RI_benchmarking_hx1")
annofpath <- file.path("eternity", "data", 
    "RI_benchmarking_resources", "gencode_v35_annotation_files",
    "gencode.v35.primary_assembly.annotation.gtf")
outdpath <- file.path("home", "metamaden", "ri_results", srrid, "superintronic")
if(!dir.exists(outdpath)){dir.create(outdpath)}

# path to features file
ftfname <- "gr-features_superintronic.rda"
ftfpath <- file.path(outdpath, ftfname)

#--------------------------
# prepare intron annotation
#--------------------------
# get features from gtf
features <- annofpath %>% collect_parts()
save(features, file = ftfpath) # save annotation file

#---------------------------
# quantify intron expression
#---------------------------
features <- get(load(ftfpath)) # load features
# iterate on chr chunks
chrv <- c(paste0("chr", seq(22)), "chrX", "chrY")
for(chr in chrv){
    message("Beginning chr chunk: ", chr, "...")
    bamfname <- paste0(srrid,'.sorted.',chr,'.bam')
    bam.fpath <- file.path(bamdpath, bamfname)
    smd <- data.frame(sample_id = srrid, bam = bam.fpath, stringsAsFactors = F)
    # get coverage from bams
    cvg <- superintronic::compute_coverage_long(smd, source = "bam")
    # get ranges overlapping gene of interest
    cvgof <- cvg %>% select(-bam) %>% join_parts(features)
    # save overlapping ranges
    cvgfname <- paste0("cov-over-ft_superintronic_",srrid,"_",chr,".rda")
    cvgfpath <- file.path(outdpath, cvgfname)
    message("Saving new results object ", cvgfpath, "...")
    save(cvgof, file = cvgfpath); message("Finished chr chunk: ", chr)}

#-------------
# save results
#-------------
# get the combined data
lf <- list.files(outdpath)
lf = lf[grepl("^cov-over-ft_superintronic.*", lf)]
tf <- do.call(rbind, lapply(lf, function(fn){
    as.data.frame(get(load(file.path(outdpath, fn))))}))
# save the combined data
tfname <- paste0("superintronic_allchr_", srrid)
rda.fpath <- file.path(outdpath, paste0(tfname, ".rda"))
tsv.fpath <- file.path(outdpath, paste0(tfname, ".tsv"))
csv.fpath <- file.path(outdpath, paste0(tfname, ".csv"))
save(tf, file = rda.fpath)
write.table(tf, sep = "\t", row.names = F, file = tsv.fpath)
write.csv(tf, row.names = F, file = csv.fpath)