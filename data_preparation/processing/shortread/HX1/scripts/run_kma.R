#!/usr/bin/env R

# Author: Sean Maden


# ENV:
# 
# 
# 1. get the dependencies from conda
# conda install -c conda-forge/label/gcc7 devtools
# conda install -c conda-forge/label/gcc7 r-reshape2=1.4.3
# conda install -c conda-forge/label/gcc7 r-dplyr=0.7.8
# conda install -c conda-forge r-data.table
# conda install -c bioconda/label/broken bowtie2
# conda install -c bioconda/label/cf201901 express
# get the python dependencies
# pip install pyfasta biopython pysam
# activate conda env
# conda activate express_kma
# 
# 2. install the lib from github, in R
# devtools::install_github("https://github.com/adamtongji/kma",
#     dependencies = FALSE)
#
# get the R lib path
# system.file("pre-process", package="kma") 
# /opt/anaconda3/envs/r351_kma/lib/R/library/kma/pre-process
# "/usr/lib64/R/library/kma/pre-process"
#
# conda install python=2.7

# 3. ENV
# setup conda env:
# conda activate py2718
# conda install -c bioconda pysam
# conda install -c anaconda biopython 
# conda install -c bioconda pyfasta


# conda activate py2718
srrid='SRR2911306'
PRE='/opt/anaconda3/envs/r351_kma/lib/R/library/kma/pre-process'
bamdpath='/eternity/data/RI_benchmarking_BAMs'
kmapreppath='eternity/data/RI_benchmarking_resources'
outdpath='/home/metamaden/ri_results/'$srrid'/kma'
sudo mkdir $outdpath

annodpath='eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/'
seqfpath=$annodpath'GRCh38.primary_assembly.genome.fa'
annofpath=$annodpath'gencode.v35.primary_assembly.annotation.gtf'

python $PRE'/generate_introns.py' --genome $seqfpath --gtf $annofpath --extend 0 --out $kmapreppath

# 4. expand and manage fastq
# srrid=SRR7410593
#find eternity/data/RI_benchmarking_fastqs -type f -iname "*SRR6026510*"
# gunzip -k eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq.gz
fqpath='eternity/data/RI_benchmarking_fastqs'
fqv=$(find $fqpath -type f -iname "*SRR2911306*")
for fn in $fqv; do sudo gzip -d $fn; done

# 5. align reads with bowtie2
# leftfastq='eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq'
# rightfastq='eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq'
# btoutpath='/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam'

annodpath='eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files'
trans_and_introns=$annodpath'/kma_reference/trans_and_introns.fa'
tipath=$annodpath'/trans_and_introns_sm'

leftfastq=$fqpath'/'$srrid'_1.fastq'
rightfastq=$fqpath'/'$srrid'_2.fastq'
btoutpath='/eternity/data/RI_benchmarking_results/kma_bams/'$srrid'_kma.bam'
sudo touch $btoutpath

fnv=$(find eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files -type f -iname "*trans_and_introns_sm*")
for fn in $fnv; do sudo chmod 777 $fn; done;
sudo chmod 777 -R eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/
sudo chmod 777 $btoutpath

sudo '/usr/bin/bowtie2' -k 200 -p 20 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x $tipath \
    -1 $leftfastq -2 $rightfastq | samtools view -Sb - > $btoutpath
# output
# 24463210 reads; of these:
#   24463210 (100.00%) were paired; of these:
#     3441596 (14.07%) aligned concordantly 0 times
#     4934735 (20.17%) aligned concordantly exactly 1 time
#     16086879 (65.76%) aligned concordantly >1 times
#     ----
#     3441596 pairs aligned concordantly 0 times; of these:
#       57275 (1.66%) aligned discordantly 1 time
#     ----
#     3384321 pairs aligned 0 times concordantly or discordantly; of these:
#       6768642 mates make up the pairs; of these:
#         5178121 (76.50%) aligned 0 times
#         482550 (7.13%) aligned exactly 1 time
#         1107971 (16.37%) aligned >1 times
# 89.42% overall alignment rate

#----------------------
# make the kma bam file
#----------------------
sudo bash -c 'cat ./kma_reference/gencode.v35.transcripts.fa ./sm/introns.fa > trans_and_introns_sm.fa'
sudo bowtie2-build --offrate 1 trans_and_introns_sm.fa trans_and_introns_sm

#------------
# run express
#------------
# download
# wget https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz

# manage paths
fapath='/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_sm.fa'
savedir='/eternity/data/RI_benchmarking_results/express_output_'$srrid'/'
sudo mkdir $savedir # make savedir

#  run express
home/metamaden/express-1.5.1-linux_x86_64/express -o $savedir $fapath $btoutpath
# returns:
# WARNING: Could not connect to update server to verify current version. Please check at the eXpress website (http://bio.math.berkeley.edu/eXpress).
# 
# 2021-Nov-12 11:21:36 - Attempting to read '/eternity/data/RI_benchmarking_hx1/SRR2911306_kma.bam' in BAM format...
# 2021-Nov-12 11:21:38 - Parsing BAM header...
# 2021-Nov-12 11:21:39 - Loading target sequences and measuring bias background...
# 2021-Nov-12 11:27:04 - Initialized 438216 targets.
# 2021-Nov-12 11:27:04 - Processing input fragment alignments...
# 2021-Nov-12 11:27:16 - Synchronized auxiliary parameter tables.

#-------------
# run kma in R
#-------------
# sudo /lib64/R/bin/$
library(kma)
srrid <- "SRR2911306"
save.dpath <- file.path("home", "metamaden", "ri_results", srrid)
resources.dpath <- file.path("eternity", "data", "RI_benchmarking_results")
xprs.dpath <- file.path(resources.dpath, paste0("express_output_",srrid))
xprs.fnv <- list.files(xprs.dpath)
xprs.results.fpath <- file.path(xprs.dpath, "results.xprs")
xprs <- read_express(xprs.results.fpath, srrid, "hx1")

# save express robject
xp.fname <- paste0("express_counts_",srrid,".rda")
xp.fpath <- file.path(save.dpath, xp.fname)
save(xprs, file = xp.fpath)

# load intron_to_trans
it.fpath <- file.path("eternity", "data", "RI_benchmarking_resources", 
    "gencode_v35_annotation_files", "sm", "intron_to_transcripts.txt")
intron_to_trans <- data.table::fread(it.fpath, data.table = FALSE)

# make new intron retention object
ir.fname <- paste0("ir-results-kma_",srrid,".rda")
ir.fpath <- file.path(save.dpath, ir.fname)
# filters
min.tpm <- 1
min.frag <- 3
# get filtered data
tpmv <- xprs$tpm; #tpmv <- tpmv[tpmv >= min.tpm]
fragv <- xprs$uniq_counts; #fragv <- fragv[fragv >= min.frag]
filt.all <- which(tpmv[,1] >= min.tpm & fragv[,1] >= min.frag)
tpmvf <- tpmv[filt.all,]
fragvf <- fragv[filt.all,]
# get intret object
ir <- newIntronRetention(targetExpression = tpmvf, unique_counts = fragvf,
    intronToUnion = intron_to_trans, groups = xprs$condition)
save(ir, file = ir.fpath)