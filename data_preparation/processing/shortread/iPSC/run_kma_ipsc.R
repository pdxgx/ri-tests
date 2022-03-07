
# 1. from command line
# get the dependencies from conda
conda install -c conda-forge/label/gcc7 devtools
conda install -c conda-forge/label/gcc7 r-reshape2=1.4.3
conda install -c conda-forge/label/gcc7 r-dplyr=0.7.8
conda install -c conda-forge r-data.table
conda install -c bioconda/label/broken bowtie2
conda install -c bioconda/label/cf201901 express
# get the python dependencies
pip install pyfasta biopython pysam

# 2. from an R session
# install the lib from github
devtools::install_github("https://github.com/adamtongji/kma",
    dependencies = FALSE)
# get the R lib path
system.file("pre-process", package="kma") 
# /opt/anaconda3/envs/r351_kma/lib/R/library/kma/pre-process

# 3. from command line
# get the intron coordinates
PRE=/opt/anaconda3/envs/r351_kma/lib/R/library/kma/pre-process
out_dir=/home/metamaden/ri_results/kma_ips/
loadpath=eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/
seqfa=echo $($loadpath'GRCh38.primary_assembly.genome.fa')
trans.gtf=$loadpath'gencode.v35.primary_assembly.annotation.gtf'
python $PRE/generate_introns.py --genome seq.fa --gtf trans.gtf --extend N --out out_dir

# 4. expand and manage fastq
# srrid=SRR7410593
#find eternity/data/RI_benchmarking_fastqs -type f -iname "*SRR6026510*"
# gunzip -k eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq.gz
fqv=$(find eternity/data/RI_benchmarking_fastqs -type f -iname "*SRR6026510*")
for fn in $fqv: do; gunzip -k $fn; done

# 5. align reads with bowtie2
leftfastq='eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq'
rightfastq='eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq'
btoutpath='/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam'
trans_and_introns='eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/kma_reference/trans_and_introns.fa'

tipath='eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_sm'
touch $btoutpath

fnv=$(find eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files -type f -iname "*trans_and_introns_sm*")
for fn in $fnv; do sudo chmod 777 $fn; done;
sudo chmod 777 -R eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/
sudo chmod 777 /eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam

sudo bowtie2 -k 200 -p 20 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x $tipath \
    -1 $leftfastq -2 $rightfastq | samtools view -Sb - > $btoutpath
# output
# 91330785 reads; of these:
#   91330785 (100.00%) were paired; of these:
#     41470026 (45.41%) aligned concordantly 0 times
#     12388302 (13.56%) aligned concordantly exactly 1 time
#     37472457 (41.03%) aligned concordantly >1 times
#     ----
#     41470026 pairs aligned concordantly 0 times; of these:
#       2019708 (4.87%) aligned discordantly 1 time
#     ----
#     39450318 pairs aligned 0 times concordantly or discordantly; of these:
#       78900636 mates make up the pairs; of these:
#         60729762 (76.97%) aligned 0 times
#         3619695 (4.59%) aligned exactly 1 time
#         14551179 (18.44%) aligned >1 times
# 66.75% overall alignment rate

#------------
# run express
#------------
# express website: https://pachterlab.github.io/eXpress/overview.html#
# download express
wget https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz
# home/metamaden/express-1.5.1-linux_x86_64/express --help

# prep the .fa file for run
# 1. combine the introns and transcripts `fa` files
sudo bash -c 'cat ./kma_reference/gencode.v35.transcripts.fa ./sm/introns.fa > trans_and_introns_sm.fa'
# 2. run alignments -- prepare bt2 file and new sample bam files....
sudo bowtie2-build --offrate 1 trans_and_introns_sm.fa trans_and_introns_sm

# manage paths
fapath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_sm.fa # trans_and_introns.fa
bamfpath=/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam
savedir=/eternity/data/RI_benchmarking_results/express_output_SRR6026510/
# make savedir
mkdir $savedir
#  run express
# express -o $outdpath/$s$outfname $fapath $bamdpath/$s$bamfname; 

# run #1 -- fails to produce results
# uses: 
# > fapath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns.fa
home/metamaden/express-1.5.1-linux_x86_64/express -o $savedir $fapath $bamfpath
# returns:
# > WARNING: Target 'chrX:49336322-49338885' exists in MultiFASTA but not alignment (SAM/BAM) file.
# then finally
# > SEVERE: Sequence for target 'chr10:95414837-95415493' not found in MultiFASTA file '/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns.fa'

# run #2
# uses:
# > fapath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns_sm.fa
home/metamaden/express-1.5.1-linux_x86_64/express -o $savedir $fapath $bamfpath
# returns:
# > 2021-Nov-09 13:54:00 - Attempting to read '/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam' in BAM format...
# > 2021-Nov-09 13:54:02 - Parsing BAM header...
# > 2021-Nov-09 13:54:03 - Loading target sequences and measuring bias background...
# > 2021-Nov-09 13:59:25 - Initialized 438216 targets.
# > 2021-Nov-09 13:59:26 - Processing input fragment alignments...
# > 2021-Nov-09 13:59:39 - Synchronized auxiliary parameter tables.
# > 2021-Nov-09 14:02:00 - Fragments Processed (/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam): 1000000        Number of Bundles: 294391.
# > 2021-Nov-09 14:02:00 - WARNING: The observed alignments appear disporportionately in the reverse-forward order (5246713 vs. 896131). If your library is strand-specific, you should use the --rf-stranded option to avoid incorrect results.
# > 2021-Nov-09 14:04:34 - Fragments Processed (/eternity/data/RI_benchmarking_results/kma_bams/SRR6026510_kma.bam): 2000000        Number of Bundles: 273363.
# ...

#-------------
# run kma in R
#-------------
# sudo /lib64/R/bin/$
library(kma)

base_dir <- system.file("example", package="kma")

xprs.dpath <- file.path("eternity", "data", "RI_benchmarking_results",
    "express_output_SRR6026510")
xprs.fnv <- list.files(xprs.dpath)
xprs.results.fpath <- file.path(xprs.dpath, "results.xprs")
xprs <- read_express(xprs.results.fpath, "SRR6026510", "ipsc")
# returns
# > Reading all eXpress data...
# > Sorting...
# > Extracting column:  tpm
# > Extracting column:  uniq_counts
names(xprs)
# [1] "tpm"         "uniq_counts" "all_data"    "condition"   "sample"
# save express robject
save.fname <- "express_counts_SRR6026510.rda"
save.dpath <- file.path("home", "metamaden", "ri_results")
save.fpath <- file.path(save.dpath, save.fname)
save(xprs, file = save.fpath)


# read in "intron_to_transcripts.txt"
i2t.fpath <- file.path("eternity", "data", "RI_benchmarking_resources", 
    "gencode_v35_annotation_files", "kma_reference", "intron_to_transcripts.txt")
intron_to_trans <- data.table::fread(i2t.fpath, data.table = FALSE)
head(intron_to_trans)
# returns:
#                     intron          target_id               gene
# 1  chr19:40733364-40736985  ENST00000597003.1  ENSG00000086544.3
# 2  chr19:40733364-40736985  ENST00000263370.3  ENSG00000086544.3
# 3   chr4:10598054-10598675  ENST00000507719.1 ENSG00000109684.15
# 4   chr4:10598054-10598675 ENST00000226951.11 ENSG00000109684.15
# 5 chr3:143796903-143832018 ENST00000316549.11 ENSG00000181804.15
# 6 chr3:143796903-143832018  ENST00000474727.2 ENSG00000181804.15
#          intron_extension strand
# 1  chr19:40733354-40736995      +
# 2  chr19:40733354-40736995      +
# 3   chr4:10598044-10598685      -
# 4   chr4:10598044-10598685      -
# 5 chr3:143796893-143832028      -
# 6 chr3:143796893-143832028      -

# make new intron retention object
ir <- newIntronRetention(xprs$tpm, intron_to_trans, xprs$condition, xprs$uniq_counts)
# returns:
# 'melting' unique counts
# computing denominator
# computing numerator
# computing retention
# 'melting' expression
# joining unique_counts and retention data
# sorting and grouping by (intron, condition)
# Warning messages:
# 1: `summarise_each_()` was deprecated in dplyr 0.7.0.
# Please use `across()` instead.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_warnings()` to see where this warning was generated.
# 2: `funs()` was deprecated in dplyr 0.8.0.
# Please use a list of either functions or lambdas:
#
#  # Simple named list:
#  list(mean = mean, median = median)
#
#  # Auto named with `tibble::lst()`:
#  tibble::lst(mean, median)
#
#  # Using lambdas
#  list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
#This warning is displayed once every 8 hours.
#Call `lifecycle::last_warnings()` to see where this warning was generated.

print(ir)
# IntronRetention object (0 introns)
# ----------------------------------------------
# Samples:        SRR6026510
# Conditions:     ipsc




