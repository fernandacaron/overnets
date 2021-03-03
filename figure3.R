rm(list=ls())

source("functions.R")

library(MASS)
library(viridis)

#This code will create plots of the expected average number of overlapping 
#species, and of the expected relationship between rho and range size overlap, 
#and between the variance of range size and variance of range size overlap.

#This function generates range overlap distribution across a set of predefined 
#meanlog values for a given sdlog, calculates the average number of overlapping 
#species by counting how many species an individual overlaps with, and taking 
#the average across all species in a matrix. The results are returned in the 
#form of points plotted in an already called plot.
meanNumOverSppFun<-function(sdlog, S, truncate, col=col){
	meanlog<-seq(from = -4, to = 0.5, by = 0.05)
	meanNumOverSpp<-matrix(ncol=2,nrow=length(meanlog))
	
	for(i in 1:length(meanlog)){
	   dat<-simRangeOver(meanlog=meanlog[i], sdlog=sdlog, S=S, truncate=truncate, 
	                     Dmin=0, Dmax=exp(1)) 
     overlap<-dat$overlap
     diag(overlap)<-0
     overlap[overlap>0]<-1
     meanNumOverSpp[i,2]<-mean(rowSums(overlap))
     meanNumOverSpp[i,1]<-dat$rho
	}
	meanNumOverSpp<-meanNumOverSpp[order(meanNumOverSpp[,1]),]
	points(meanNumOverSpp[,1], meanNumOverSpp[,2], col=col, pch=16)
}


#This function plots points corresponding to the expected mean range overlap 
#for a rho value if sdlog is given, or the expected variance of range overlap 
#for a value of the variance of range size if meanlog is given.
plotRho<-function(meanlog, sdlog, S, truncate, col=col){
  if(missing(meanlog)){
    meanlog<-seq(from = -4, to = 0.5, by = 0.012)
    rhoMean<-matrix(ncol=2, nrow=length(meanlog))
    
    for(i in 1:length(meanlog)){
      dat<-simRangeOver(meanlog=meanlog[i], sdlog=sdlog, S=S, truncate=truncate, 
                        Dmin=0, Dmax=exp(1)) 
      overlap<-dat$overlap
      range.sizes<-dat$range.sizes
      rhoMean[i,1]<-dat$rho
      diag(overlap)<-NA
      overlap<-overlap[overlap > 0]
      rhoMean[i,2]<-mean(log(overlap), na.rm = TRUE)
    }
    points(rhoMean[,1], rhoMean[,2], col=col, pch = 16, cex=0.7)
    } else {
    sdlog<-seq(from = 0, to = 10, by = 0.02)
    rhoSd<-matrix(ncol=2, nrow=length(sdlog))
    
    for(i in 1:length(sdlog)){
      dat<-simRangeOver(meanlog=meanlog, sdlog=sdlog[i], S=S, truncate=truncate, 
                        Dmin=0, Dmax=exp(1)) 
      overlap<-dat$overlap
      range.sizes<-dat$range.sizes
      rhoSd[i,1]<-var(log(range.sizes))
      diag(overlap)<-NA
      overlap<-overlap[overlap > 0]
      rhoSd[i,2]<-var(as.vector(log(overlap)), na.rm = TRUE)
    }
    points(rhoSd[,1], rhoSd[,2], col=col, pch = 16, cex = 0.7)
  }
}

S<-10000

#FigA
pdf("Figure3.pdf")
layout(matrix(1:9, ncol=3, byrow=TRUE), widths = c(1, 1, 0.4))
par(mar = c(4, 4, 4, 1))
plot(-1, ylim=c(0, S), xlim=c(-4, 0.5), type="n", xlab=expression(rho), 
     ylab="Average number of overlapping species", frame.plot=F)

rbPal <- colorRampPalette(c(viridis(10)), alpha=TRUE)
col <- rbPal(10)[as.numeric(cut(1:10, breaks = 10))]

meanNumOverSppFun(sdlog=0.1, S=S, truncate=TRUE, col=col[1])
meanNumOverSppFun(sdlog=0.5, S=S, truncate=TRUE, col=col[2])
meanNumOverSppFun(sdlog=0.7, S=S, truncate=TRUE, col=col[3])
meanNumOverSppFun(sdlog=1, S=S, truncate=TRUE, col=col[4])
meanNumOverSppFun(sdlog=1.5, S=S, truncate=TRUE, col=col[5])
meanNumOverSppFun(sdlog=1.7, S=S, truncate=TRUE, col=col[6])
meanNumOverSppFun(sdlog=2, S=S, truncate=TRUE, col=col[7])
meanNumOverSppFun(sdlog=2.5, S=S, truncate=TRUE, col=col[8])
meanNumOverSppFun(sdlog=2.7, S=S, truncate=TRUE, col=col[9])
meanNumOverSppFun(sdlog=3, S=S, truncate=TRUE, col=col[10])

#FigB
plot(-1, ylim=c(0, S), xlim=c(-4, 0.5), type="n", xlab=expression(rho), 
     ylab="", frame.plot=F)

rbPal <- colorRampPalette(c(viridis(10)), alpha=TRUE)
col <- rbPal(10)[as.numeric(cut(1:10, breaks = 10))]

meanNumOverSppFun(sdlog=0.1, S=S, truncate=FALSE, col=col[1])
meanNumOverSppFun(sdlog=0.5, S=S, truncate=FALSE, col=col[2])
meanNumOverSppFun(sdlog=0.7, S=S, truncate=FALSE, col=col[3])
meanNumOverSppFun(sdlog=1, S=S, truncate=FALSE, col=col[4])
meanNumOverSppFun(sdlog=1.5, S=S, truncate=FALSE, col=col[5])
meanNumOverSppFun(sdlog=1.7, S=S, truncate=FALSE, col=col[6])
meanNumOverSppFun(sdlog=2, S=S, truncate=FALSE, col=col[7])
meanNumOverSppFun(sdlog=2.5, S=S, truncate=FALSE, col=col[8])
meanNumOverSppFun(sdlog=2.7, S=S, truncate=FALSE, col=col[9])
meanNumOverSppFun(sdlog=3, S=S, truncate=FALSE, col=col[10])

par(mar = c(4, 3, 4, 2))

color.bar <- function(lut, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = "") {
  scale = (length(lut) )/(max - min)
  
  plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = title)
  axis(2, ticks, las = 1)
  for (i in 1:(length(lut) )) {
    y = (i - 1)/scale + min
    rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
  }
}
color.bar(col, min = 0, max = 3, title = "")
mtext("sdlog", cex=0.8)


S<-1000

#FigC
par(mar = c(4, 4, 4, 1))

plot(-1, ylim=c(-5,0), xlim=c(-4, 0.5), type="n", xlab=expression(rho), 
     ylab="Mean range overlap", frame.plot=F)

plotRho(sdlog=0.1, S=S, truncate=TRUE, col=adjustcolor(viridis(20)[6], alpha.f = 0.8))
plotRho(sdlog=0.8, S=S, truncate=TRUE, col=adjustcolor(viridis(20)[10], alpha.f = 0.8))
plotRho(sdlog=1.5, S=S, truncate=TRUE, col=adjustcolor(viridis(20)[12], alpha.f = 0.8))
plotRho(sdlog=2.3, S=S, truncate=TRUE, col=adjustcolor(viridis(20)[14], alpha.f = 0.8))
plotRho(sdlog=3, S=S, truncate=TRUE, col=adjustcolor(viridis(20)[16], alpha.f = 0.8))

#FigD
plot(-1, ylim=c(-5,0), xlim=c(-4, 0.5), type="n", xlab=expression(rho), 
     ylab="", frame.plot=F)

plotRho(sdlog=0.1, S=S, truncate=FALSE, col=adjustcolor(viridis(20)[6], alpha.f = 0.8))
plotRho(sdlog=0.8, S=S, truncate=FALSE, col=adjustcolor(viridis(20)[10], alpha.f = 0.8))
plotRho(sdlog=1.5, S=S, truncate=FALSE, col=adjustcolor(viridis(20)[12], alpha.f = 0.8))
plotRho(sdlog=2.3, S=S, truncate=FALSE, col=adjustcolor(viridis(20)[14], alpha.f = 0.8))
plotRho(sdlog=3, S=S, truncate=FALSE, col=adjustcolor(viridis(20)[16], alpha.f = 0.8))

par(mar = c(4, 3, 4, 2))

color.bar(viridis(20)[6:16], min = 0.1, max = 3, title = "")
mtext("sdlog", cex=0.8)

#FigE
par(mar = c(4, 4, 4, 1))

plot(-1, ylim=c(0,90), xlim=c(0, 115), type="n", xlab="Variance range size", 
     ylab="Variance range overlap", frame.plot=F)

plotRho(meanlog=-4, S=S, truncate=TRUE, col=adjustcolor(plasma(20)[6], alpha.f = 0.8))
plotRho(meanlog=-2.9, S=S, truncate=TRUE, col=adjustcolor(plasma(20)[10], alpha.f = 0.8))
plotRho(meanlog=-1.75, S=S, truncate=TRUE, col=adjustcolor(plasma(20)[12], alpha.f = 0.8))
plotRho(meanlog=-0.6, S=S, truncate=TRUE, col=adjustcolor(plasma(20)[14], alpha.f = 0.8))
plotRho(meanlog=0.5, S=S, truncate=TRUE, col=adjustcolor(plasma(20)[16], alpha.f = 0.8))
abline(b=1, a=0)

#FigF
plot(-1, ylim=c(0,90), xlim=c(0, 115), type="n", xlab="Variance range size", 
     ylab="", frame.plot=F)

plotRho(meanlog=-4, S=S, truncate=FALSE, col=adjustcolor(plasma(20)[6], alpha.f = 0.8))
plotRho(meanlog=-2.9, S=S, truncate=FALSE, col=adjustcolor(plasma(20)[10], alpha.f = 0.8))
plotRho(meanlog=-1.75, S=S, truncate=FALSE, col=adjustcolor(plasma(20)[12], alpha.f = 0.8))
plotRho(meanlog=-0.6, S=S, truncate=FALSE, col=adjustcolor(plasma(20)[14], alpha.f = 0.8))
plotRho(meanlog=0.5, S=S, truncate=FALSE, col=adjustcolor(plasma(20)[16], alpha.f = 0.8))
abline(b=1, a=0)

par(mar = c(4, 3, 4, 2))

color.bar(plasma(20)[6:16], min = -4, max = 0.5, title = "")
mtext("meanlog", cex=0.8)

dev.off()

