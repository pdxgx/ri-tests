#!/usr/bin/env R

# Author: Sean Maden
#
# Correlate intron features with SR IR tool expressions.

library(ggplot2)
library(ggpubr)
library(gridExtra)

srrid <- "SRR6026510"
run.handle <- "ipsc"
plot.title <- "iPSC"

#-----------
# load data
#-----------
# tsv called RIs
tsv.fname <- "called_RI_data_summary_iPSCfeatureannotated.tsv"
tsv.bi <- read.table(tsv.fname, sep = "\t", header = T)

#--------------
# get plot data
#--------------
# make dfp
cnv.comp1 <- c("width", "intron_position_in_tx", "total_overlapping_features", 
               "max_features_per_base", "X._bases_overlapped")
tsvf <- tsv.bi
cnv.comp2 <- colnames(tsvf)[grepl(".*weighted_median$|lwm$", colnames(tsvf))]
compv <- paste0(cnv.comp2, ":", rep(cnv.comp1, each = length(cnv.comp2)))
dfcor <- do.call(rbind, lapply(seq(length(compv)), function(ii){
  var1 <- gsub(":.*", "", compv[ii]); var2 <- gsub(".*:", "", compv[ii])
  tsvff <- tsvf[tsvf[,var2] > 0,]; cti <- cor.test(tsvff[,var1], tsvff[,var2])
  data.frame(cor.pval = format(cti$p.value, scientific = T), 
             cor.rho = round(cti$estimate, 3), 
             var1 = var1, var2 = var2)}))

# format vars
dfp <- dfcor
dfp[grepl("interest", dfp$var1),]$var1 <- "IntEREst"
dfp[grepl("iread", dfp$var1),]$var1 <- "iREAD"
dfp[grepl("irfinders", dfp$var1),]$var1 <- "IRFinder-S"
dfp[grepl("kma", dfp$var1),]$var1 <- "KMA"
dfp[grepl("superintronic", dfp$var1),]$var1 <- "superintronic"
dfp[grepl("position_in_tx", dfp$var2),]$var2 <- "Transcript position"
dfp[grepl("max_features_per_base", dfp$var2),]$var2 <- "Max feat. per base"
dfp[grepl("total_overlapping", dfp$var2),]$var2 <- "Total overlapping feat."
dfp[grepl("width", dfp$var2),]$var2 <- "Intron length (bp)"
dfp[grepl("X._bases", dfp$var2),]$var2 <- "Bases overlapping feat."
dfp$`Rho` <- dfp$cor.rho
dfp$label <- as.character(round(as.numeric(dfp$cor.rho), 2))

#----------------------------------------------------
# make new heatmap figure -- all and filtered introns
#----------------------------------------------------
# get heatmap objects
ggcor <- ggplot(dfp, aes(x = var1, y = var2, fill = `Rho`)) + 
  geom_tile(color = "navyblue", size = 1) + theme_bw() + ylab("Intron property") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.1),
        axis.title.x = element_blank(), legend.position = "none") +
  geom_text(aes(label = label), color = "#6B6B6B") + ggtitle(plot.title) + 
  scale_fill_gradient2(low = "#520B71", mid = "white", midpoint = 0, high = "#E68400")

plot.legend <- ggplot(dfp, aes(x = var1, y = var2, fill = Rho)) + geom_tile(color = "white") + 
  scale_fill_gradient2(low = "#520B71", mid = "white", midpoint = 0, high = "#E68400",
                       name = expression(rho)) +
  theme(legend.title=element_text(size=14, face = "italic"))
plot.legend <- get_legend(plot.legend)

# save new figs
plot.fname <- paste0("mcor_intron-features_sr-ir-5tools_", srrid,"-",run.handle)
lm <- matrix(c(rep(1, 6), 2),nrow = 1)
ylab.str <- paste0(paste0(rep(" ", 15), collapse = ""), "Intron feature")
# save new pdf
pdf(paste0(plot.fname, ".pdf"), 6, 4)
grid.arrange(ggcor, as_ggplot(plot.legend), nrow = 1, layout_matrix = lm)
dev.off()
# save new png
png(paste0(plot.fname, ".png"), width = 5,  height = 3.1, res = 500, units = "in")
grid.arrange(ggcor, as_ggplot(plot.legend), nrow = 1, layout_matrix = lm)
dev.off()