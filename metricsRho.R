rm(list=ls())

library(sna)
library(viridis)

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

D<-1
rho<-seq(from = 0.001, to = 1, by = 0.002)
meanlog<-rho*D
S<-1000
nsim<-100

calc_TRUE<-function(meanlog){
  met1<-list()
  met1$dg<-met1$bt<-met1$cl<-matrix(nrow=S)
  dg<-bt<-cl<-matrix(nrow=S, ncol=nsim)
  for(j in 1:nsim){
  overlap<-simRangeOver(meanlog=meanlog, S=S, truncate=TRUE)
  diag(overlap)<-NA
  
  dg[,j]<-calcMetrics(overlap)$degree
  bt[,j]<-calcMetrics(overlap)$betweeness
  cl[,j]<-calcMetrics(overlap)$closeness
  }
  met1$dg<-sort(rowMeans(dg), decreasing=T)
  met1$bt<-sort(rowMeans(bt), decreasing=T)
  met1$cl<-sort(rowMeans(cl), decreasing=T)
  return(met1)
}


calc_FALSE<-function(meanlog){
  met2<-list()
  met2$dg<-met2$bt<-met2$cl<-matrix(nrow=S)
  dg<-bt<-cl<-matrix(nrow=S, ncol=nsim)
  for(j in 1:nsim){
    overlap<-simRangeOver(meanlog=meanlog, S=S, truncate=FALSE)
    diag(overlap)<-NA
    
    dg[,j]<-calcMetrics(overlap)$degree
    bt[,j]<-calcMetrics(overlap)$betweeness
    cl[,j]<-calcMetrics(overlap)$closeness
  }
  met2$dg<-sort(rowMeans(dg), decreasing=T)
  met2$bt<-sort(rowMeans(bt), decreasing=T)
  met2$cl<-sort(rowMeans(cl), decreasing=T)
  return(met2)
}


met1<-list()
met1$dg<-met1$bt<-met1$cl<-matrix(nrow=S, ncol=length(meanlog))
for (j in 1:length(meanlog)){
  input<-list()
  for(i in 1:nsim){
    input[[i]]<-c(meanlog[j])
  }
  xx<-lapply(input, calc_TRUE)
  dg<-bt<-cl<-matrix(nrow=S, ncol=nsim)
  for(i in 1:nsim){
    dg[,i]<-xx[[i]]$dg
    bt[,i]<-xx[[i]]$bt
    cl[,i]<-xx[[i]]$cl
  }
  met1$dg[,j]<-sort(rowMeans(dg), decreasing=T)
  met1$bt[,j]<-sort(rowMeans(bt), decreasing=T)
  met1$cl[,j]<-sort(rowMeans(cl), decreasing=T)
}


met2<-list()
met2$dg<-met2$bt<-met2$cl<-matrix(nrow=S, ncol=length(meanlog))
for (j in 1:length(meanlog)){
  input<-list()
  for(i in 1:nsim){
    input[[i]]<-c(meanlog[j])
  }
  xx<-lapply(input, calc_TRUE)
  dg<-bt<-cl<-matrix(nrow=S, ncol=nsim)
  for(i in 1:nsim){
    dg[,i]<-xx[[i]]$dg
    bt[,i]<-xx[[i]]$bt
    cl[,i]<-xx[[i]]$cl
  }
  met2$dg[,j]<-sort(rowMeans(dg), decreasing=T)
  met2$bt[,j]<-sort(rowMeans(bt), decreasing=T)
  met2$cl[,j]<-sort(rowMeans(cl), decreasing=T)
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
