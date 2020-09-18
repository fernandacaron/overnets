rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

final1<-final2<-as.data.frame(matrix(ncol = 3, nrow = 100))
colnames(final1)<-colnames(final2)<-c("rho", "mean.range.size", "mean.range.overlap")

rho<-seq(from = 0.1, to = 1, by = 0.1)
meanlog<-1/rho
S<-1000
for (i in 1:100){
  x<-simRangeOver(meanlog=meanlog[i], S=S, truncate=TRUE)
  overlap<-x
  
  #range size stats
  mean.range.size<-mean(diag(x))
  
  #range overlap stats
  diag(x)<-NA
  #x<-matrix(x[x > 0], ncol = 1) #pq só os maiores que 0? não teria que considerar
                                #todos?
  mean.range.overlap <- mean(x, na.rm = TRUE)

  res<-cbind(rho[i], mean.range.size, mean.range.overlap)
  
  final1[i,]<-res
}
for (i in 1:100){
  x<-simRangeOver(meanlog=meanlog[i], S=S, truncate=FALSE)
  overlap<-x
  
  #range size stats
  mean.range.size<-mean(diag(x))

    #range overlap stats
  diag(x) <- NA
  #x <- matrix(x[x > 0], ncol = 1)
  mean.range.overlap <- mean(x, na.rm = TRUE)

  res<-cbind(rho[i], mean.range.size, mean.range.overlap)
  
  final2[i,]<-res
}

layout(matrix(1:3, ncol=2))

plot(mean.range.size ~rho, data=final1, xlab = expression(rho), main="Range size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(mean.range.overlap ~rho, data=final1, xlab = expression(rho), main="Range overlap size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(mean.range.overlap ~ mean.range.size, data=final1, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap ~ mean.range.size, data=final1))

