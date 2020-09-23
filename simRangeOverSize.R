rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

rho<-seq(from = 0.01, to = 1, by = 0.01)
meanlog<-rho #rho = meanlog/D(=1)
S<-1000

final1<-final2<-as.data.frame(matrix(ncol = 3, nrow = length(rho)))
colnames(final1)<-colnames(final2)<-c("rho", "mean.range.size", "mean.range.overlap")

for (i in 1:length(rho)){
  overlap<-simRangeOver(meanlog=meanlog[i], S=S, truncate=TRUE)
  
  #range size stats
  mean.range.size<-mean(diag(overlap))
  
  #range overlap stats
  diag(overlap)<-NA
  mean.range.overlap <- mean(overlap, na.rm = TRUE)

  res<-cbind(rho[i], mean.range.size, mean.range.overlap)
  
  final1[i,]<-res
}

for (i in 1:length(rho)){
  overlap<-simRangeOver(meanlog=meanlog[i], S=S, truncate=FALSE)
  
  #range size stats
  mean.range.size<-mean(diag(overlap))

  #range overlap stats
  diag(overlap) <- NA
  mean.range.overlap <- mean(overlap, na.rm = TRUE)

  res<-cbind(rho[i], mean.range.size, mean.range.overlap)
  
  final2[i,]<-res
}

layout(matrix(1:2, ncol=2))

plot(mean.range.overlap~mean.range.size, data=final1, ylim=c(0.3,0.8),  xlab="Range size", main="Truncated", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final1))
abline(a=0, b=1, col="red")

plot(mean.range.overlap~mean.range.size, data=final2, xlim=c(1,5), xlab="Range size", main="", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final2))
abline(a=0, b=1, col="red")

plot(mean.range.size~rho, data=final1, xlab = expression(rho), main="Range size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.size~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(mean.range.overlap~rho, data=final1, xlab = expression(rho), main="Range overlap size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.overlap~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)


