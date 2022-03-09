#!/usr/bin/env R

# Author: Sean Maden
# 
# Format the superintronic results. Bind and summarize adjacent regions.
#

# set the run identifiers
srrid <- "SRR2911306"
run.handle <- "hx1"

#----------
# load data
#----------
si.fname <- paste0("superintronic_allchr_",srrid, ".rda")
si <- get(load(si.fname))

#-------------------------
# filter duplicate regions
#-------------------------
# get the feature labels/id
cname <- "score"
si$intronid <- paste0(si$seqnames, ":", si$start, "-", si$end)

# get duplicated features
dup <- which(duplicated(si$intronid))
dup.id <- unique(si$intronid[dup])
# get new df for duplicated ids
which.notdup <- si$intronid %in% dup.id & !duplicated(si$intronid)
df.notdup <- si[which.notdup,]
# retain max observed expr by dup id
dfply <- si[si$intronid %in% dup.id,]
eval.str <- paste0("dfply <- dfply %>% group_by(intronid) ",
                   "%>% summarise(",cname," = max(",cname,", na.rm = T))")
eval(parse(text = eval.str)); 
dfply <- dfply[order(match(dfply$intronid, df.notdup$intronid)),]
df.notdup[,cname] <- dfply[,2]
# bind results
sif <- rbind(si[!si$intronid %in% dup.id,], df.notdup)

#---------------------
# get granges and save
#---------------------
si.gr <- makeGRangesFromDataFrame(sif, keep.extra.columns = T)

# save granges
si.gr.fname <- paste0("granges_superintronic_", srrid,"-",run.handle,".rda")
save(si.gr, file = si.gr.fname)

# save df
si.gr.df <- as.data.frame(si.gr, stringsAsFactors = FALSE)
si.gr.df.fname <- paste0("df-granges_superintronic_", srrid,"-",run.handle, ".csv")
write.csv(si.gr.df, file = si.gr.df.fname)