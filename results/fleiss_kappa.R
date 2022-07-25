library("irr")
HX1 <- read.table("called_RI_data_summary_HX1featureannotated_GCcontent.tsv",sep="\t",header=TRUE)
introns <- unique(HX1[,"intron"])
HX1.dat <- data.frame(matrix(nr=length(introns),nc=8))
for(i in 1:length(introns)) {
	HX1.dat[i,] <- apply(HX1[which(HX1[,"intron"]==introns[i]),3:10],2,max)
}
kappam.fleiss(HX1.dat>0)
# Fleiss' Kappa for m Raters
#
# Subjects = 12687 
#   Raters = 8 
#    Kappa = 0.145 
#
#        z = 86.5 
#  p-value = 0 
iPSC <- read.table("called_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv",sep="\t",header=TRUE)
introns <- unique(iPSC[,"intron"])
iPSC.dat <- data.frame(matrix(nr=length(introns),nc=8))
for(i in 1:length(introns)) {
	iPSC.dat[i,] <- apply(iPSC[which(iPSC[,"intron"]==introns[i]),3:10],2,max)
}
kappam.fleiss(iPSC.dat>0)
# Fleiss' Kappa for m Raters
#
# Subjects = 7900 
#   Raters = 8 
#    Kappa = 0.068 
#
#        z = 32 
#  p-value = 0 
combined <- rbind(iPSC.dat, HX1.dat)
kappam.fleiss(combined>0)
# Fleiss' Kappa for m Raters
#
# Subjects = 20587 
#   Raters = 8 
#    Kappa = 0.113 
#
#        z = 85.6 
#  p-value = 0 
