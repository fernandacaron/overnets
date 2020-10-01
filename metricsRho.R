rm(list=ls())

library(sna)

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

D<-1
rho<-seq(from = 0.1, to = 1, by = 0.1)
meanlog<-rho*D
S<-1000

met1<-list()
met1$dg<-met1$bt<-met1$cl<-matrix(nrow=S, ncol=length(rho))
for (i in 1:length(rho)){
  nsim<-10
  dg<-bt<-cl<-matrix(nrow=S, ncol=nsim)
  for(j in 1:nsim){
    overlap<-simRangeOver(meanlog=meanlog[1], S=S, truncate=TRUE)
    diag(overlap)<-NA
    
    dg[,j]<-calcMetrics(overlap)$degree
    bt[,j]<-calcMetrics(overlap)$betweeness
    cl[,j]<-calcMetrics(overlap)$closeness
  }
  met1$dg[,i]<-sort(rowMeans(dg), decreasing=T)
  met1$bt[,i]<-sort(rowMeans(bt), decreasing=T)
  met1$cl[,i]<-sort(rowMeans(cl), decreasing=T)
}

layout(matrix(1:4, ncol=4, byrow=TRUE), widths = c(1, 1, 1, 0.4))
par(mar = c(4, 4, 4, 1))
col <- magma(10)[as.numeric(cut(rho, breaks = 10))]

plot(met1$dg[,1], ylab="Degree", xlab="rank", ylim=range(met1$dg), type="n") 
for (i in 1:length(rho)){
  lines(met1$dg[,i], col=col[i])
}

plot(met1$bt[,1], ylab="Betweenness", xlab="rank", ylim=range(met1$bt), type="n") 
for (i in 1:100){
  lines(met1$bt[,i], col=col[i])
}

plot(met1$cl[,1], ylab="Closeness", xlab="rank", ylim=range(met1$cl), type="n") 
for (i in 1:100){
  lines(met1$cl[,i], col=col[i])
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

color.bar(col, min = 0, max = 1, title = "")
mtext(expression(rho))
