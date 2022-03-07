#!/usr/bin/env bash

# Author: Sean Maden
# Run the iREAD software
# 
# NOTES:
# 
# * get gff annotation file:
#       -- https://www.gencodegenes.org/human/
#       -- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz
# ENV:
# 
# activate the conda env with python 2###
# > conda activate py2718_iread

# manage variables and paths
srrid=SRR2911306
bamdpath='/eternity/data/RI_benchmarking_hx1'
bedpath='/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.introns.bed' 
ireadpath=iread.py #/home/metamaden/iread/iread.py
outputpath=$srrid'_iread'
sudo mkdir $outputpath

#----------
# run iread
#----------

# iterate on chr chunks, implicit sudo access
cd iread
sudo chmod -R 777 .
for chri in $chrv
    do
        echo 'Beginning chr '$chri
        bamfpath=$bamdpath/$srrid'.sorted.'$chri'.nochr.bam'
        numreads=$(samtools view -c -F 260 $bamfpath)
        python iread.py $bamfpath $bedpath -o $outputpath -t $numreads
        echo 'Finished chr '$chri
    done

#--------------------------
# load iread results into R
#--------------------------
sudo /lib64/R/bin/R

srrid <- "SRR2911306"
newtname <- paste0(srrid,"_allchr_iread.txt")

readdpath <- "SRR2911306_iread"
lf <- list.files(readdpath)
lf <- lf[grepl(paste0(srrid, "\\.sorted\\.chr.*"), lf)]
lf <- lf[!grepl(".*_intron_reads.*", lf)]


tf <- do.call(rbind, lapply(lf, function(fn){
    message(fn); read.table(fn, sep = "\t", header = T)}))
write.table(tf, file = newtname, row.names = F)








