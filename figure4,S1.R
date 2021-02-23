rm(list=ls())

library(sna)
library(viridis)

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

netStats <- function(fixed.sdlog = TRUE, S, truncate = TRUE) {
	nStats <- function(meanlog, sdlog, S, Dmin, Dmax, truncate = TRUE) {
		nreps <- 100
		res <- list()
		bt <- cl <- dg  <- matrix(ncol = S, nrow = nreps)
		for (i in 1:nreps) {
			x <- simRangeOver(meanlog=meanlog, sdlog=sdlog, S=S, Dmin=Dmin, 
								Dmax=Dmax, truncate = truncate)$overlap
			diag(x)<-NA
			x[x > 0] <- 1
			bt[i, ] <- sort(betweenness(x, gmode = "graph"), decreasing = TRUE)
			cl[i, ] <- sort(closeness(x, gmode = "graph"), decreasing = TRUE)
			dg[i, ] <- sort(degree(x, gmode = "graph"), decreasing = TRUE)
		}
		res$bt <- colMeans(bt)
		res$cl <- colMeans(cl)
		res$dg <- colMeans(dg)
		res
	}
	nreps <- 100
	meanlog <- seq(from = -4, to = 0.5, by = 0.045)
	sdlog<-seq(from = 0, to = 10, by = 0.1)
	Dmin <- 0
	Dmax <- exp(1)
	if(fixed.sdlog==TRUE){
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:nreps) {
	print(i)
		x <- nStats(meanlog[i], sdlog=1, S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	} else {
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:nreps) {
	print(i)
		x <- nStats(meanlog=0, sdlog=sdlog[i], S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	}
}

meanlogTr<-netStats(fixed.sdlog=TRUE, S=200, truncate=TRUE)
meanlogN<-netStats(fixed.sdlog=TRUE, S=200, truncate=FALSE)
sdlogTr<-netStats(fixed.sdlog=FALSE, S=200, truncate=TRUE)
sdlogN<-netStats(fixed.sdlog=FALSE, S=200, truncate=FALSE)

pdf("figures/Figure4.pdf", width=9, height=7)
layout(matrix(1:8, ncol=4, byrow=T), widths = c(1, 1, 1, 0.4))
meanlog <- seq(from = -4, to = 0.5, by = 0.045)
sdlog<-seq(from = 0, to = 10, by = 0.1)
rbPal1 <- colorRampPalette(c(plasma(20)[6:16]), alpha=TRUE)
col1 <- rbPal1(101)[as.numeric(cut(1:length(meanlog), breaks = 101))]

rbPal2 <- colorRampPalette(c(viridis(20)[6:16]), alpha=TRUE)
col2 <- rbPal2(101)[as.numeric(cut(1:length(sdlog), breaks = 101))]

par(mar=c(4, 4, 2, 2))

plot(meanlogN$dg[1,], ylab="Degree", xlab="", ylim=range(meanlogN$dg), type="n") 
for (i in 1:length(meanlog)){
  lines(meanlogN$dg[i,], col=col1[i])
}

plot(meanlogN$bt[1,], ylab="Betwenness", xlab="", ylim=range(meanlogN$bt), type="n") 
for (i in 1:length(meanlog)){
  lines(meanlogN$bt[i,], col=col1[i])
}

plot(meanlogN$cl[1,], ylab="Closeness", xlab="", ylim=range(meanlogN$cl), type="n") 
for (i in 1:length(rho_mean_NT)){
  lines(meanlogN$cl[i,], col=col1[i])
}

color.bar <- function(lut, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = "") {
  scale = (length(lut) - 1)/(max - min)
  
  #   dev.new(width=1.75, height=5)\n
  plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = title)
  axis(2, ticks, las = 1)
  for (i in 1:(length(lut) - 1)) {
    y = (i - 1)/scale + min
    rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
  }
}
par(mar = c(3, 3, 3, 2))

y1<-round(meanlog, digits=2)
color.bar(col1, min = min(y1), max = max(y1), title = "")
mtext("meanlog", cex=0.8)

plot(sdlogN$dg[1,], ylab="Degree", xlab="rank", ylim=range(sdlogN$dg), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogN$dg[i,], col=col2[i])
}

plot(sdlogN$bt[1,], ylab="Betwenness", xlab="rank", ylim=range(sdlogN$bt), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogN$bt[i,], col=col2[i])
}

plot(sdlogN$cl[1,], ylab="Closeness", xlab="rank", ylim=range(sdlogN$cl), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogN$cl[i,], col=col2[i])
}

par(mar = c(3, 3, 3, 2))

y2<-round(sdlog, digits=1)
color.bar(col2, min = min(y2), max = max(y2), title = "")
mtext("sdlog", cex=0.8)

dev.off()
 

pdf("figures/FigureS1.pdf", width=9, height=7)
layout(matrix(1:8, ncol=4, byrow=T), widths = c(1, 1, 1, 0.4))
meanlog <- seq(from = -4, to = 0.5, by = 0.045)
sdlog<-seq(from = 0, to = 10, by = 0.1)
rbPal1 <- colorRampPalette(c(plasma(20)[6:16]), alpha=TRUE)
col1 <- rbPal1(101)[as.numeric(cut(1:length(meanlog), breaks = 101))]

rbPal2 <- colorRampPalette(c(viridis(20)[6:16]), alpha=TRUE)
col2 <- rbPal2(101)[as.numeric(cut(1:length(sdlog), breaks = 101))]

par(mar=c(4, 4, 2, 2))

plot(meanlogTr$dg[1,], ylab="Degree", xlab="", ylim=range(meanlogTr$dg), type="n") 
for (i in 1:length(meanlog)){
  lines(meanlogTr$dg[i,], col=col1[i])
}

plot(meanlogTr$bt[1,], ylab="Betwenness", xlab="", ylim=range(meanlogTr$bt), type="n") 
for (i in 1:length(meanlog)){
  lines(meanlogTr$bt[i,], col=col1[i])
}

plot(meanlogTr$cl[1,], ylab="Closeness", xlab="", ylim=range(meanlogTr$cl), type="n") 
for (i in 1:length(rho_mean_NT)){
  lines(meanlogTr$cl[i,], col=col1[i])
}

par(mar = c(3, 3, 3, 2))

y1<-round(meanlog, digits=2)
color.bar(col1, min = min(y1), max = max(y1), title = "")
mtext("meanlog", cex=0.8)

plot(sdlogTr$dg[1,], ylab="Degree", xlab="rank", ylim=range(sdlogTr$dg), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogTr$dg[i,], col=col2[i])
}

plot(sdlogTr$bt[1,], ylab="Betwenness", xlab="rank", ylim=range(sdlogTr$bt), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogTr$bt[i,], col=col2[i])
}

plot(sdlogTr$cl[1,], ylab="Closeness", xlab="rank", ylim=range(sdlogTr$cl), type="n") 
for (i in 1:length(sdlog)){
  lines(sdlogTr$cl[i,], col=col2[i])
}

par(mar = c(3, 3, 3, 2))

y2<-round(sdlog, digits=1)
color.bar(col2, min = min(y2), max = max(y2), title = "")
mtext("sdlog", cex=0.8)

dev.off()