genes.fname.hx1 <- "HX1_validatedONLY_RI_genes_07-16-2022_07.46.06.tsv"
genes.fname.ipsc <- "iPSC_validatedONLY_RI_genes_07-16-2022_07.46.06.tsv"
lgenes <- list()
lgenes[["HX1"]] <- read.table(genes.fname.hx1, sep = "\t", header = T)
lgenes[["iPSC"]] <- read.table(genes.fname.ipsc, sep = "\t", header = T)
toolv <- c("IntEREst", "iREAD", "IRFinder.S", "superintronic", "kma", "rMATS", "MAJIQ", "SUPPA2")
cnv <- colnames(lgenes[[1]]); 
dfp <- do.call(rbind, lapply(names(lgenes), function(samplei){
  mgenei <- lgenes[[samplei]]
  do.call(rbind, lapply(toolv, function(tooli){
    cnvi <- c(cnv[grepl(tooli, cnv)], "intron", "gene_name")
    mf <- mgenei[,cnvi]; intronv <- unique(mf$intron)
    do.call(rbind, lapply(intronv, function(introni){
      mfi <- mf[mf$intron==introni,]; ni <- names(mfi)
      # message(introni)
      if(length(mfi) > 0){
        tmi <- ifelse(mfi[grepl("true_positives", ni)] == 1, "TP",
                      ifelse(mfi[grepl("false_positives", ni)] == 1, "FP", 
                             ifelse(mfi[grepl("false_negatives", ni)] == 1, "FN", "TN")))
        data.frame("tool" = tooli, "intron" = introni, "sample" = samplei,
                   "tmetric" = as.character(tmi), "gene" = mfi$gene_name)
      }
    }))
    #message(tooli)
  }))
  #message(samplei)
}))
#genev <- unique(dfp$gene)
genev <- c("AP1G2","LBR","SRSF7","IGSF8","FAHD2B")
lgg <- lapply(genev, function(genei){
  dfpi <- dfp[dfp$gene==genei,]
  samplev <- unique(dfpi$sample)
  dfpi$start <- as.numeric(gsub(".*:|-.*", "", dfpi$intron))
  dfpi$end <- as.numeric(gsub(".*-", "", dfpi$intron))
  intronv <- unique(dfpi$intron)
  startv <- as.numeric(gsub(".*:|-.*", "", intronv))
  intronv <- intronv[order(startv)]
  dfpi$intron.lab <- 0
  for(ii in seq(length(intronv))){
    dfpi[dfpi$intron==intronv[ii],]$intron.lab <- ii
  }
  dfpi <- unique(dfpi)
  retmat <- matrix(nrow=length(toolv)*length(lgenes), ncol=max(dfpi$intron.lab), dimnames=list(rep(toolv,length(lgenes)), 1:max(dfpi$intron.lab)))
  for (i in 1:max(dfpi$intron.lab)) {
  	print(paste(genei, dfpi[which(dfpi$intron.lab == i), "tmetric"]))
  	for (j in 1:length(toolv)) {
  		for(k in 1:length(lgenes)) {
  			res.ijk <- unique(dfpi[which(dfpi$tool == toolv[j] & dfpi$intron.lab == i & dfpi$sample == names(lgenes)[k]), "tmetric"])
  			if (length(res.ijk) > 0) {
		  		retmat[j+length(toolv)*(k-1),i] <- res.ijk
		  	}
 		}
  	}
  }
  return(retmat) })
  

plot_mat <- function(x, y, mat=NULL, label="", split=9, col.highlight=0,size=NULL) {
	if (is.null(mat)) {
		return()
	}
	x.max <- par("usr")[2]
	x.min <- par("usr")[1]
	x.width <- (x.max-x.min)/75
	y.min <- par("usr")[3]
	y.max <- par("usr")[4]
	y.width <- (y.max-y.min)/75
	x.dim <- dim(mat)[1]
	y.dim <- dim(mat)[2]
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
			adj <- j %/% split
			x1 <- x+i*x.width
			x2 <- x1+x.width
			y1 <- y-(1+j+adj)*y.width
			y2 <- y-(2+j+adj)*y.width
			
			if (is.na(mat[j,i])) {
				rect(x1, y1, x2, y2, col="white", border=bord)
			}
			else if (mat[j,i] == "FP") {
				rect(x1, y1, x2, y2, col="#d8788a", border=bord)
			}
			else if (mat[j,i] == "TP") {
				rect(x1, y1, x2, y2, col="#9db92c", border=bord)
			}
			else if (mat[j,i] == "TN") {
				rect(x1, y1, x2, y2, col="#ffc18f", border=bord)
			}
			else if (mat[j,i] == "FN") {
				rect(x1, y1, x2, y2, col="#6d9cc6", border=bord)
			}
		}
	}
	text(x+x.width,y-y.width, label, cex=0.75, adj=0, font=3)		
	for (j in 1:x.dim) {
		adj <- j %/% split
		yj <- y-(1.5+j+adj)*y.width
		text(x+x.width/2, yj, j %% split + j %/% split, cex=0.5, adj=0.5)
	}
	text(x, y-0.5*split*y.width, "HX1", cex=0.75, srt=90, adj=1, pos=2)
	text(x, y-1.5*split*y.width, "iPSC", cex=0.75, srt=90, adj=1, pos=2)
}

par(mar=c(0.1,0.1,0.1,0.1))
plot(0:1,0:1,type="n",axes=FALSE, xlab="", ylab="")
plot_mat(0,1,mat=lgg[[1]], label=genev[1], col.highlight=17) # AP1G2
plot_mat(0.46,1,mat=lgg[[2]], label=genev[2], col.highlight=5) # LBR
plot_mat(0.06,0.69,mat=lgg[[3]], label=genev[3], col.highlight=7) # SRSF7
plot_mat(0.3,0.69,mat=lgg[[4]], label=genev[4], col.highlight=5) # IGSF8
plot_mat(0.47,0.69,mat=lgg[[5]], label=genev[5], col.highlight=1) # FAHD2B
legend(0.01,0.38,legend=c("TP", "FN", "FP", "TN", "No coverage"), fill=c("#9db92c","#6d9cc6","#d8788a","#ffc18f","white"),border="lightgray", ncol=4, cex=0.75)
