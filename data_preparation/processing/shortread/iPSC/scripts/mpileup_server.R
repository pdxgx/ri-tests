#!/usr/bin/env R

# Author: Sean Maden
#
# Get median short-read RNA-seq coverage by gene using samtools mpileup.
# This result informs the coverage filters prior to analyses of retained
# introns.
#

#library(rtracklayer)
library(data.table)

# set the run identifiers
srrid <- "SRR6026510"
run.handle <- "ipsc"

# set annotation version
gencode_version <- "35"

#-----------
# load data
#-----------
# define paths
results_dpath <- "/home/metamaden/ri_results"
bamfname = paste0(srrid, '.sorted.bam')
# write gene info in bed format
bed.fpath <- "/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3.bed"
bed <- read.table(bed.fpath)
# filter bed file
bed.genes <- bed[bed[,8]=="gene",]
bed.genes.fname <- "gencode_genes.v35.annotation.bed"
bed.genes.fpath <- file.path(results_dpath, bed.genes.fname)
write.table(bed.genes, file = bed.genes.fpath, sep = " ", 
            row.names = F, col.names = F, quote = FALSE)

#----------------------------------------------
# get gene pileup summaries -- short read data
#----------------------------------------------
# define vars
bam.dpath <- file.path("eternity", "data", "RI_benchmarking_BAMs")
bam.fpath <- file.path(bam.dpath, paste0(srrid, '.sorted.bam')) # bam read path
mpiletmp.fpath <- file.path(results_dpath, paste0('mpiletemp_',srrid,'.tsv')) # mpiletmp path
bedtmp.fpath <- file.path(results_dpath, paste0('bedtmp_',srrid,'.bed')) # bedtmp path
num.round <- 2 # value round num digits
cnv <- c("gene.id", "chr", "start", "end", "length", "mean", 
         "median", "sd", "q0", "q25", "q50", "q75", "q100")
mp.fname <- paste0("mpile-sstat_gencode-v",gencode_version,
    "_",srrid,"-",run.handle, ".csv")
mp.fpath <- file.path(results_dpath, mp.fname) # results table write path
fwrite(matrix(cnv, nrow = 1), file = mp.fpath, sep = ",",
       append = F, col.names = F, row.names = F) # save new results table

# read in genes bedfile
bed.genes <- fread(file = bed.genes.fpath, sep = " ", data.table = F) 
uchrv <- unique(bed.genes$V1); uchrv <- as.character(uchrv[grepl("chr", uchrv)])
samtools.path <- "/usr/local/bin/samtools"; numgenes.chunk <- 2000
t1 <- Sys.time()
for(ii in seq(length(uchrv))){
  chri <- uchrv[ii]; message("Beginning chr: ", chri)
  bedf.all <- bed.genes[bed.genes[,1] == chri,];itot <- nrow(bedf.all)
  for(jj in seq(1, itot, numgenes.chunk)){
    end.index <- ifelse(jj+numgenes.chunk-1>itot,itot,jj+numgenes.chunk-1)
    bedf.chunk <- bedf.all[seq(jj, end.index),]
    write.table(as.matrix(bedf.chunk), file = bedtmp.fpath, sep = " ", 
                row.names = F, col.names = F, quote = F) # bedtmp
    sys.str <- paste0(c(samtools.path, "mpileup", "-r", chri, "-l", bedtmp.fpath,
                        bam.fpath, ">", mpiletmp.fpath), collapse = " ")
    system(sys.str) # run samtools, get coverage/bp for the gene
    mpu.chunk <- read.table(file = mpiletmp.fpath, sep = "\t")
    # parse results
    for(ri in seq(nrow(bedf.chunk))){
      message("Writing gene results...")
      line <- bedf.chunk[ri,]; geneidi <- line[4]
      coordv <- as.numeric(c(line[2], line[3]))
      coorddiff <- abs(coordv[1]-coordv[2]) # get gene length
      # get summary stats
      filt.mpu <- mpu.chunk[,2] >= min(coordv) & mpu.chunk[,2] <= max(coordv)
      mpu.ri <- mpu.chunk[filt.mpu,]; puv <- as.numeric(mpu.ri[,4])
      meani <- round(mean(puv, na.rm = T), num.round)
      mediani <- round(median(puv, na.rm = T), num.round)
      sdi <- round(sd(puv, na.rm = T), num.round)
      quantilei <- unlist(lapply(quantile(puv, seq(0,1,0.25)), 
                                 round, num.round))
      mdat <- matrix(c(as.character(geneidi[1,1]), chri, coordv[1], 
                       coordv[2], coorddiff, meani, mediani, sdi, 
                       quantilei), nrow = 1)
      fwrite(mdat, file = mp.fpath, append = T, col.names = F, 
             row.names = F, sep = ",") # write new data
    };message("Finished chunk ",jj,", time: ", Sys.time() - t1)
  };message("Finished chr: ", ii, ", time: ", Sys.time() - t1)
}

#---------------------------
# pileup summaries --introns
#---------------------------
# read in intron bed file
ibed.fname <- "introns.bed"; ibed <- read.table(ibed.fname)
ibed[,1] <- paste0("chr", ibed[,1]); dim(ibed)
# define vars
bamfname <- paste0(srrid, '.sorted.bam'); num.round <- 2
mpiletmp <- 'mpiletemp.tsv'; bedtmp <- 'bedtmp.bed' 
# make new return table
# table name
mpfname <- paste0("mpile-sstat_gencode-v",gencode_version,
                  "-introns_", srrid, "-", run.handle, ".tsv")
cnv <- c("chr", "start", "end", "strand", "gene.id", "mean", 
         "median", "sd", "q0", "q25", "q50", "q75", "q100")
# save new table colnames
fwrite(matrix(cnv, nrow = 1), file = mpfname, sep = ",", 
       append = F, col.names = F, row.names = F)
# begin process
t1 <- Sys.time(); uchrv <- unique(bed.genes$V1); numgenes.chunk <- 1000
for(ii in seq(length(uchrv))){
  chri <- uchrv[ii]; message("Beginning chr: ", chri)
  bedf.all <- bed.genes[bed.genes[,1] == chri,];itot <- nrow(bedf.all)
  for(jj in seq(1, itot, numgenes.chunk)){
    end.index <- jj+numgenes.chunk-1; end.index <- ifelse(end.index>itot,itot,end.index)
    rv <- seq(jj, end.index); bedf.chunk <- bedf.all[rv,] # get row indices
    write.table(as.matrix(bedf.all[rv,]), file = bedtmp, sep = " ", 
                row.names = F, col.names = F, quote = F) # bedtmp
    sys.str <- paste0(c("samtools", "mpileup", "-r", chri, "-l", bedtmp, 
                        bamfname, ">", mpiletmp), collapse = " ")
    system(sys.str); mpu.chunk <- read.table(file = mpiletmp, sep = "\t") # get pileup data
    for(ri in seq(nrow(bedf.chunk))){ # parse results
      line <- bedf.chunk[ri,]; geneidi <- line[4]; coordv <- as.numeric(c(line[2], line[3]))
      ibed.genef <- ibed[ibed$V1==chri & ibed[,2] >= min(coordv) & ibed[,3] <= max(coordv),]
      message("Parsing ",nrow(ibed.genef), " introns for gene ", geneidi, "...")
      if(nrow(ibed.genef) > 0){message("Writing intron results...")
        icoordv <- as.numeric(unlist(c(ibed.genef[2], ibed.genef[3])))
        puv <- as.numeric(mpu.chunk[mpu.chunk[,2] >= min(icoordv) & mpu.chunk[,2] <= max(icoordv),4])
        meani <- round(mean(puv, na.rm = T), num.round)
        mediani <- round(median(puv, na.rm = T), num.round)
        sdi <- round(sd(puv, na.rm = T), num.round)
        quantilei <- unlist(lapply(quantile(puv, seq(0,1,0.25)), round, num.round))
        mdat <- matrix(c(chri, line[2], line[3], line[6], 
                         geneidi, meani, mediani, sdi, quantilei), nrow = 1)
        fwrite(mdat, file = mpfname, append = T, col.names = F, 
               row.names = F, sep = ",")}
    }; message("Finished chunk ",jj,", time: ", Sys.time() - t1)
  }; message("Finished chr: ", ii, ", time: ", Sys.time() - t1)
}

#----------------------------------------------
# get gene pileup summaries -- iterate on genes
#----------------------------------------------
# define vars
mpiletmp <- 'mpiletemp.tsv'
bedtmp <- 'bedtmp.bed'
num.round <- 2
# make new return table
cnv <- c("gene.id", "chr", "start", "end", "length", "mean", 
                    "median", "sd", "q0", "q25", "q50", "q75", "q100")
mpfname <- paste0("mpile-sstat_gencode-v38_", srrid, "-", run.handle, ".tsv")
fwrite(matrix(cnv, nrow = 1), file = mpfname, sep = ",",
       append = F, col.names = F, row.names = F)
# read in genes bedfile
bed.genes <- fread(file = bed.genes.fname, sep = " ", data.table = F) 
t1 <- Sys.time()
for(ii in seq(nrow(bed.genes))){
  message("Beginning gene ", ii); line <- bed.genes[ii,]
  # write bedtmp
  write.table(matrix(line, nrow = 1), file = bedtmp, sep = " ", 
              row.names = F, col.names = F, quote = F)
  # run samtools, get coverage/bp for the gene
  sys.str <- paste0(c("samtools", "mpileup", "-r", chri, "-l", bedtmp, 
                      bamfname, ">", mpiletmp), collapse = " ")
  system(sys.str)
  # parse results
  chri = line[1]; geneidi = line[4] # get gene point info
  coordstart = as.numeric(line[2]); coordend = as.numeric(line[3])
  coorddiff = abs(coordstart-coordend) # get gene length
  # read in results
  file.conn <- file(mpiletmp); mpu <- read.table(file = mpiletmp, sep = "\t")
  # get summary stats
  puv <- as.numeric(mpu[,4]); meani <- round(mean(puv, na.rm = T), num.round)
  mediani <- round(median(puv, na.rm = T), num.round); sdi <- round(sd(puv, na.rm = T), num.round)
  quantilei <- unlist(lapply(quantile(puv, seq(0,1,0.25)), round, num.round))
  # return data
  mdat <- matrix(c(geneidi, chri, coordstart, coordend, coorddiff, 
                   meani, mediani, sdi, quantilei), nrow = 1)
  fwrite(mdat, file = mpfname, append = T, col.names = F, 
         row.names = F, sep = ",")
  message("Finished gene ", ii, ", time: ", Sys.time() - t1)
}

#--------------------
# analysis of results
#--------------------
dat <- fread(file = mpfname, sep = ",", header = T, data.table = F)
boxplot(as.numeric(dat$median))
plot(density(as.numeric(dat$median)))