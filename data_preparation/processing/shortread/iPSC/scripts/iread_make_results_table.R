# sudo /lib64/R/bin/R

# Author: Sean Maden
#
# Load results from iREAD into R. 
#

# get vars
srrid <- "SRR2911306"
newtname <- paste0(srrid,"_allchr_iread.txt")
readdpath <- "SRR2911306_iread"
lf <- list.files(readdpath)
lf <- lf[grepl(paste0(srrid, "\\.sorted\\.chr.*"), lf)]
lf <- lf[!grepl(".*_intron_reads.*", lf)]

# load data
tf <- do.call(rbind, lapply(lf, function(fn){
    message(fn); read.table(fn, sep = "\t", header = T)}))

# save new table
write.table(tf, file = newtname, row.names = F)