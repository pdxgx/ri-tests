#!/usr/bin/evn R

# Author: Sean Maden
#
# Run IntEREst to quantify intron expression from chr-chunked 
# BAM file. This uses the GFF3 reference file for GENCODE v35.
# 
# Run R on the remote server with sudo access:
# > sudo /lib64/R/bin/R

library(IntEREst)

# set run id
srrid="SRR2911306"

#----------
# load data
#----------
# define paths
results.dpath <- file.path("home", "metamaden", "ri_results", srrid)
interest.dpath <- file.path(results.dpath, "interest")
if(!dir.exists(results.dpath)){dir.create(results.dpath)}
if(!dir.exists(interest.dpath)){dir.create(interest.dpath)}

# get the bam file chunks path
bam.dpath <- file.path("eternity", "data", "RI_benchmarking_hx1")
lfv <- list.files(bam.dpath)
which.chr.bam <- grepl("\\.chr.*", lfv) & grepl("\\.bam$", lfv) &
    !grepl(".*\\.nochr.*", lfv)
lfv <- lfv[which.chr.bam]

# load the reference object
ref.fpath <- file.path("home", "metamaden", "ri_results", 
    "interest_ipsc", "gen35gff3_interest-ref.rda")
ref <- get(load(ref.fpath))

#---------------
# process chunks
#---------------
# process each bam file chunk
for(fn in lfv){
    message("Working on file ", fn)
    bam.fpath <- file.path(bam.dpath, fn)
    # get the output file paths
    chri <- gsub(".*sorted\\.|\\.bam", "", bam.fpath)
    out.fname <- paste0("interest_",chri,"_",srrid,".tsv")
    log.fname <- paste0("log_",chri,"_",srrid,".txt")
    out.fpath <- file.path(interest.dpath, out.fname)
    log.fpath <- file.path(interest.dpath, log.fname)
    message("Running interest...")
    results <- IntEREst::interest(
        bamFileYieldSize=10000, junctionReadsOnly=FALSE,
        bamFile=bam.fpath, isPaired=TRUE, isPairedDuplicate=FALSE,
        isSingleReadDuplicate=NA, reference=ref, 
        referenceGeneNames=ref[,"collapsed_gene_id"],
        referenceIntronExon=ref[,"int_ex"],
        outFile=out.fpath, logFile=log.fpath, method="IntRet", clusterNo=1, 
        returnObj=TRUE, scaleLength= TRUE, scaleFragment= TRUE)
    message("Finished with file ", fn)
}

#------------------
# get results table
#------------------
# read results chunks into single table file
# interest.dpath <- "."; srrid="SRR6026510"
fnv <- list.files(interest.dpath)
fnv <- fnv[grepl("\\.tsv$", fnv) & grepl("^interest_chr.*", fnv)]
tf <- do.call(rbind, lapply(fnv, function(fn){
    chri <- gsub("(^interest_|_SRR.*)", "", fn)
    tsvi <- data.table::fread(file.path(interest.dpath, fn), 
        sep="\t", data.table = F)
    tsvif <- tsvi[tsvi[,1] == chri,]; return(tsvif)}))
dim(tf)
tf$intronid <- paste0(tf$chr, ":", tf$begin, "-", tf$end)

# save results table
tfname <- paste0("interest_allchr_", srrid)
rda.fpath <- file.path(interest.dpath, paste0(tfname, ".rda"))
txt.fpath <- file.path(interest.dpath, paste0(tfname, ".txt"))
csv.fpath <- file.path(interest.dpath, paste0(tfname, ".csv"))
save(tf,file = rda.fpath)
write.csv(tf, file = csv.fpath)
write.table(tf, file = txt.fpath)
