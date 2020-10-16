rm(list=ls())

library(sna)
library(viridis)
library(svMisc)

setwd("~/overnets")

source("functions.R")

meanlog<-seq(from = -5, to = 1, by = 0.06)
S<-100
reps<-100

met1<-list()
rho1<-numeric()
met1$dg<-met1$bt<-met1$cl<-matrix(nrow=S, ncol=length(meanlog))
for (j in 1:length(meanlog)){
  progress(j, max.value=length(meanlog))
  dg<-bt<-cl<-matrix(nrow=S, ncol=reps)
  for(i in 1:reps){
    overlap<-simRangeOver(meanlog=meanlog[j], S=S, truncate=TRUE)$over
    rho1[j]<-mean(diag(overlap))
    diag(overlap)<-NA
    #overlap[overlap > 0] <- 1
    
    dg[,i]<-sort(calcMetrics(overlap)$degree, decreasing=T)
    bt[,i]<-sort(calcMetrics(overlap)$betweeness, decreasing=T)
    cl[,i]<-sort(calcMetrics(overlap)$closeness, decreasing=T)
  }
  met1$dg[,j]<-rowMeans(dg)
  met1$bt[,j]<-rowMeans(bt)
  met1$cl[,j]<-rowMeans(cl)
}

rho2<-numeric()
met2<-list()
met2$dg<-met2$bt<-met2$cl<-matrix(nrow=S, ncol=length(meanlog))
for (j in 1:length(meanlog)){
  progress(j, max.value=length(meanlog))
  dg<-bt<-cl<-matrix(nrow=S, ncol=reps)
  for(i in 1:reps){
    overlap<-simRangeOver(meanlog=meanlog[j], S=S, truncate=FALSE)$over
    rho2[j]<-mean(diag(overlap))
    diag(overlap)<-NA
    #overlap[overlap > 0] <- 1
    
    dg[,i]<-sort(calcMetrics(overlap)$degree, decreasing=T)
    bt[,i]<-sort(calcMetrics(overlap)$betweeness, decreasing=T)
    cl[,i]<-sort(calcMetrics(overlap)$closeness, decreasing=T)
  }
  met2$dg[,j]<-rowMeans(dg)
  met2$bt[,j]<-rowMeans(bt)
  met2$cl[,j]<-rowMeans(cl)
}


pdf("figures/Figure4.pdf", width=8, height=4)
layout(matrix(1:4, ncol=4, byrow=TRUE), widths = c(1, 1, 1, 0.4))
rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col1 <- rbPal1(101)[as.numeric(cut(1:length(rho1), breaks = 101))]

plot(met1$dg[,1], ylab="Degree", xlab="rank", ylim=range(met1$dg), type="n") 
for (i in 1:length(rho1)){
  lines(met1$dg[,i], col=col1[i])
}

plot(met1$bt[,1], ylab="Betweenness", xlab="rank", ylim=range(met1$bt), type="n") 
for (i in 1:length(rho1)){
  lines(met1$bt[,i], col=col1[i])
}

plot(met1$cl[,1], ylab="Closeness", xlab="rank", ylim=range(met1$cl), type="n") 
for (i in 1:length(rho1)){
  lines(met1$cl[,i], col=col1[i])
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

y1<-round(rho1, digits=2)
color.bar(col1, min = min(y1), max = max(y1), title = "")
mtext(expression(rho))
dev.off()


pdf("figures/FigureS2.pdf", width=8, height=4)
layout(matrix(1:4, ncol=4, byrow=TRUE), widths = c(1, 1, 1, 0.4))
rbPal2 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col2 <- rbPal2(101)[as.numeric(cut(1:length(rho2), breaks = 101))]

plot(met2$dg[,1], ylab="Degree", xlab="rank", ylim=range(met2$dg), type="n") 
for (i in 1:length(rho2)){
  lines(met2$dg[,i], col=col2[i])
}

plot(met2$bt[,1], ylab="Betweenness", xlab="rank", ylim=range(met2$bt), type="n") 
for (i in 1:length(rho2)){
  lines(met2$bt[,i], col=col2[i])
}

plot(met2$cl[,1], ylab="Closeness", xlab="rank", ylim=range(met2$cl), type="n") 
for (i in 1:length(rho2)){
  lines(met2$cl[,i], col=col2[i])
}

y2<-round(rho2, digits=2)
color.bar(col2, min = min(y2), max = max(y2), title = "")
mtext(expression(rho))
dev.off()