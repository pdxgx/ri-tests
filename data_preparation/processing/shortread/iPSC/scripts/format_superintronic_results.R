#!/usr/bin/env R

# Author: Sean Maden
# Format the superintronic results. 
# Bind and summarize adjacent regions.

srrid <- "SRR6026510"; run.handle <- "ipsc"

#----------
# load data
#----------
si.fname <- paste0("superintronic_allchr_",srrid, ".rda")
si <- get(load(si.fname))
dim(si) 
# [1] 47282793       12

#-------------------------
# filter duplicate regions
#-------------------------
cname <- "score"

# get the feature labels/id
si$intronid <- paste0(si$seqnames, ":", si$start, "-", si$end)

# get duplicated features
dup <- which(duplicated(si$intronid))
dup.id <- unique(si$intronid[dup])
length(dup.id) 
# [1] 6838923

# get new df for duplicated ids
which.notdup <- si$intronid %in% dup.id & !duplicated(si$intronid)
df.notdup <- si[which.notdup,]
dim(df.notdup) 
# [1] 6838923      13
# retain max observed expr by dup id
dfply <- si[si$intronid %in% dup.id,]; dim(dfply) 
# [1] 14724394       13
eval.str <- paste0("dfply <- dfply %>% group_by(intronid) ",
                   "%>% summarise(",cname," = max(",cname,", na.rm = T))")
eval(parse(text = eval.str)); 
dim(dfply) 
# [1] 6838923       2
dfply <- dfply[order(match(dfply$intronid, df.notdup$intronid)),]
identical(dfply$intronid, df.notdup$intronid) 
# TRUE
df.notdup[,cname] <- dfply[,2]
# bind results
sif <- rbind(si[!si$intronid %in% dup.id,], df.notdup)
dim(si) 
# [1] 47282793       13
dim(sif) 
# [1] 39397322       13

#---------------------
# get granges and save
#---------------------
si.gr <- makeGRangesFromDataFrame(sif, keep.extra.columns = T)

# save granges
si.gr.fname <- paste0("granges_superintronic_",
                      srrid,"-",run.handle,".rda")
save(si.gr, file = si.gr.fname)

# save df
si.gr.df <- as.data.frame(si.gr, stringsAsFactors = FALSE)
si.gr.df.fname <- paste0("df-granges_superintronic_",
                         srrid,"-",run.handle, ".csv")
write.csv(si.gr.df, file = si.gr.df.fname)
