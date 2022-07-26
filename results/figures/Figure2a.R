#!/usr/bin/env Rscript --vanilla

#HELPER FUNCTION AND SCRIPT TO GENERATE FIGURE 2a

# plot_mat() - helper function to overlay long-read data matrix for a given transcript onto plot
# input parameters as follows:
#   pt: numerical value specifying where along x-axis values to plot corresponding data
#   ht: numerical value specifying y-axis coordinate to extend line from plot
#   mat: matrix of 0/1/NA values describing read-level data for a transcript; each row is a read, each col is an intron
#   col.highlight: integer value specifying which col in matrix to treat as relevant intron 
#   kink: 2-value vector specifying (+x, +y) offsets to introduce kink into display line
#   label: character value with which to label input matrix in plot
#   data: stored output from hist() function
plot_mat <- function(pt, ht=NULL, mat=NULL, col.highlight=0, kink=c(0,0), label="", data) {
	p <- which.min(abs(data$mids - pt))
	x <- data$mids[p]
	y <- data$counts[p]
	if (is.null(ht)) {
		ht <- y
	}
	lines(c(x,x,x+kink[1]), c(y, ht,ht+kink[2]), lty="dotted")
	points(x, y, pch=19, col="lightgray")
	points(x, y, pch=1)
	x <- x+kink[1]
	ht <- ht+kink[2]
	if (is.null(mat)) {
		return()
	}
	x.width <- (par("usr")[2]-par("usr")[1])/75
	y.width <- (par("usr")[4]-par("usr")[3])/75
	x.dim <- dim(mat)[1]
	y.dim <- dim(mat)[2]
	col.landing <- mean(col.highlight)
	for (i in c(setdiff(1:y.dim, col.highlight), col.highlight)) {
		if (i==0) {
			next
		}
		if (i %in% col.highlight) {
			bord <- "black"
		}
		else {
			bord <- "lightgray"
		}
		for (j in 1:x.dim) {
			if (is.na(mat[j,i])) {
				rect(x+(i-col.landing-0.5)*x.width, ht+(x.dim-j)*y.width, x+(i+1-col.landing-0.5)*x.width, ht+(x.dim+1-j)*y.width, col="white", border=bord)
			}
			else if (mat[j,i] > 0) {
				rect(x+(i-col.landing-0.5)*x.width, ht+(x.dim-j)*y.width, x+(i+1-col.landing-0.5)*x.width, ht+(x.dim+1-j)*y.width, col="#d9f0a3", border=bord)
			}
			else {
				rect(x+(i-col.landing-0.5)*x.width, ht+(x.dim-j)*y.width, x+(i+1-col.landing-0.5)*x.width, ht+(x.dim+1-j)*y.width, col="#006034", border=bord)
			}
		}
	}
	text(x+(y.dim/2 - col.landing+0.5)*x.width,ht+(x.dim+1)*y.width, label, cex=0.5, adj=0.5)
}

#READ AND PROCESS DATA
#read in iPSC csv data, with columns corresponding to 'intron', 'transcript', 'ir_confidence', and 'position'
iPSC.path <- readline(prompt = "iPSC filepath: ") # processed_tx_df_iPSC_02-21-2022_21.03.50.csv
iPSC <- read.csv(as.character(iPSC.path)) 
#exclude from distribution introns that are entirely spliced out
iPSC.results <- iPSC[which(iPSC[,3] > 0), 3]
#generate 100-bin histogram of iPSC intron persistence data
data.iPSC <- hist(iPSC.results, plot=FALSE, n=100)

#GENERATE FIGURE
plot(data.iPSC$mids,data.iPSC$counts,type="l", lwd=2, ylab="Introns (#)", xlab=bquote("Persistence ("~P["i,t"]~")"), xaxt="n", cex.lab=1.4)
axis(1, at=0:5/5,labels=c(">0",1:5/5))
legend("topright", legend=c("Spliced", "Retained", "No coverage"), fill=c("#006034","#d9f0a3", "white"), border=c("black", "black", "gray"), ncol=2)

#ENST00000446856.5
mat <- matrix(0, nr=16,nc=15)
mat[,10] <- 1
plot_mat(1, ht=715, mat, col.highlight=10, kink=c(-0.1,50), label="ENST00000446856.5", data=data.iPSC)

#ENST00000539455.5
mat <- matrix(nr=7,nc=4)
mat[1,] <- c(0,0,0,1)
mat[2,] <- c(0,0,0,1)
mat[3,] <- c(0,0,0,1)
mat[4,] <- c(0,0,0,1)
mat[5,] <- c(0,0,1,1)
mat[6,] <- c(NA,NA,0,1)
mat[7,] <- c(NA,NA,0,1)
plot_mat(0.85, ht=100, mat, col.highlight=4, kink=c(0.05,50), label="ENST00000539455.5", data=data.iPSC)

#ENST00000426395.7
mat <- matrix(0, nr=12, nc=5)
mat[,4] <- 1
mat[c(9,11,12),2] <- 1
mat[10,5] <- 1
plot_mat(0.78, ht=425, mat, col.highlight=4, label="ENST00000426395.7", data=data.iPSC)

#ENST00000329908.12
mat <- matrix(nr=12,nc=9)
mat[1:2,] <- 0
mat[3,] <- c(0,0,0,0,0,0,0,1,0)
mat[4,] <- c(0,0,0,0,0,0,0,1,0)
mat[5,] <- c(0,0,0,0,0,0,0,1,0)
mat[6,] <- c(0,0,0,0,0,0,0,1,0)
mat[7,] <- c(0,0,0,0,1,0,0,1,0)
mat[8,] <- c(0,0,0,0,0,0,1,1,0)
mat[9,] <- c(0,0,0,0,0,0,1,1,0)
mat[10,] <- c(0,0,0,0,1,0,0,1,0)
mat[11,] <- c(0,0,0,0,0,0,1,1,0)
mat[12,] <- c(0,0,0,0,1,1,0,1,0)
plot_mat(0.64, ht=40, mat, col.highlight=8, kink=c(0.06,35), label="ENST00000329908.12", data=data.iPSC)

#ENST00000645141.1
mat <- matrix(0, nr=7, nc=8)
mat[,4] <- 1
mat[4:5,3] <- 1
mat[6:7,5] <- 1
mat[6,6] <- 1
mat[c(2,7),7] <- 1
mat[c(3,5:7),8] <- 1
plot_mat(0.55, ht=450, mat, col.highlight=4, label="ENST00000645141.1", data=data.iPSC)

#ENST00000580026.5
mat <- matrix(nr=12,nc=3)
mat[,1:2] <- 1
mat[,3] <- 0
mat[1,] <- 0
mat[2,2] <- 0
plot_mat(0.44, ht=40, mat, col.highlight=1, label="ENST00000580026.5", data=data.iPSC)

#ENST00000429589.5
mat <- matrix(nr=15, nc=9)
mat[1,] <- c(0,0,0,0,0,0,0,0,1)
mat[2,] <- c(0,0,0,0,0,0,0,0,1)
mat[3,] <- c(0,0,0,0,0,0,0,0,1)
mat[4,] <- c(0,0,0,0,0,1,0,0,0)
mat[5,] <- c(0,0,0,0,1,0,0,0,1)
mat[6,] <- c(NA,NA,NA,0,0,0,1,0,1)
mat[7,] <- c(NA,NA,NA,1,0,0,0,0,1)
mat[8,] <- c(NA,0,1,0,0,0,1,1,1)
mat[9,] <- c(NA,NA,NA,0,1,1,1,0,1)
mat[10,] <- c(NA,NA,NA,0,1,0,1,1,1)
mat[11,] <- c(NA,NA,NA,NA,NA,NA,NA,0,0)
mat[12,] <- c(NA,NA,NA,0,1,1,1,1,1)
mat[13,] <- c(NA,1,1,0,1,1,1,1,1)
mat[14,] <- c(NA,0,1,1,1,1,1,1,1)
mat[15,] <- c(NA,1,0,1,1,1,1,1,1)
plot_mat(0.28, ht=300, mat, col.highlight=9, kink=c(0.1,100), label="ENST00000429589.5", data=data.iPSC)

#ENST00000292807.9
mat <- matrix(nr=7,nc=11)
mat[1,] <- c(0,0,0,1,0,0,0,0,0,0,1)
mat[2,] <- c(0,0,0,1,1,0,1,0,0,0,1)
mat[3,] <- c(0,0,0,1,0,0,1,0,0,1,1)
mat[4,] <- c(0,1,0,1,0,0,1,0,0,1,1)
mat[5,] <- c(0,0,0,1,1,0,1,0,0,1,1)
mat[6,] <- c(0,0,0,1,1,0,1,0,0,1,1)
mat[7,] <- c(0,1,0,1,0,0,1,0,0,1,1)
plot_mat(0.17, ht=775, mat, col.highlight=2, kink=c(0.025,25), label="ENST00000292807.9", data=data.iPSC)

#ENST00000544500.5
mat <- matrix(nr=8,nc=8)
mat[1:2,] <- 0
mat[3,] <- c(NA,0,0,0,0,0,0,0)
mat[4,] <- c(0,0,0,0,0,1,1,0)
mat[5,] <- c(NA,1,0,0,0,1,1,0)
mat[6,] <- c(NA,0,0,0,0,1,1,0)
mat[7,] <- c(NA,NA,NA,0,0,1,1,0)
mat[8,] <- c(NA,NA,0,1,1,1,1,0)
plot_mat(0.06,ht=1050,mat,col.highlight=2, label="ENST00000544500.5", data=data.iPSC)
