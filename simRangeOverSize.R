rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

rho<-seq(from = 0.001, to = 1, by = 0.002)
meanlog<-rho #rho = meanlog/D(=1)
S<-1000

final1<-final2<-as.data.frame(matrix(ncol = 9, nrow = length(rho)))
colnames(final1)<-colnames(final2)<-c("rho", "mean.range.size", "mean.range.overlap",
                                      "variance.range.size", "variance.range.overlap",
                                      "skewness.range.size", "skewness.range.overlap",
                                      "kurtosis.range.size", "kurtosis.range.overlap")

for (i in 1:length(rho)){
  overlap<-simRangeOver(meanlog=meanlog[i], S=S, truncate=TRUE)
  
  #range size stats
  mean.range.size<-mean(diag(overlap))
  variance.range.size<-var(diag(overlap))
  skewness.range.size<-skewness(diag(overlap))
  kurtosis.range.size<-kurtosis(diag(overlap))
  
  #range overlap stats
  diag(overlap)<-NA
  mean.range.overlap<-mean(overlap, na.rm = TRUE)
  variance.range.overlap<-var(as.vector(overlap), na.rm = TRUE)
  skewness.range.overlap<-skewness(as.vector(overlap), na.rm = TRUE)
  kurtosis.range.overlap<-kurtosis(as.vector(overlap), na.rm = TRUE)
  
  res<-cbind(rho[i], mean.range.size, mean.range.overlap, variance.range.size, 
             variance.range.overlap, skewness.range.size, skewness.range.overlap,
             kurtosis.range.size, kurtosis.range.overlap)
  
  final1[i,]<-res
}

for (i in 1:length(rho)){
  overlap<-simRangeOver(meanlog=meanlog[i], S=S, truncate=FALSE)
  
  #range size stats
  mean.range.size<-mean(diag(overlap))
  variance.range.size<-var(diag(overlap))
  skewness.range.size<-skewness(diag(overlap))
  kurtosis.range.size<-kurtosis(diag(overlap))
  
  #range overlap stats
  diag(overlap) <- NA
  mean.range.overlap <- mean(overlap, na.rm = TRUE)
  variance.range.overlap<-var(as.vector(overlap), na.rm = TRUE)
  skewness.range.overlap<-skewness(as.vector(overlap), na.rm = TRUE)
  kurtosis.range.overlap<-kurtosis(as.vector(overlap), na.rm = TRUE)
  
  res<-cbind(rho[i], mean.range.size, mean.range.overlap, variance.range.size, 
             variance.range.overlap, skewness.range.size, skewness.range.overlap,
             kurtosis.range.size, kurtosis.range.overlap)
  
  final2[i,]<-res
}

#Mean and Variance
layout(matrix(1:4, ncol=2))

plot(mean.range.overlap~mean.range.size, data=final1, xlab="Range size", main="Mean", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final1))

plot(mean.range.overlap~mean.range.size, data=final2, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final2))

plot(variance.range.overlap~variance.range.size, data=final1, xlab="Range size", main="Variance", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(variance.range.overlap~variance.range.size, data=final1))

plot(variance.range.overlap~variance.range.size, data=final2, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
abline(lm(variance.range.overlap~variance.range.size, data=final2))


#Skewness and Kurtosis 
layout(matrix(1:4, ncol=2))

plot(skewness.range.overlap~skewness.range.size, data=final1, xlab="Range size", main="Skewness", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(skewness.range.overlap~skewness.range.size, data=final1))

plot(skewness.range.overlap~skewness.range.size, data=final2, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
abline(lm(skewness.range.overlap~skewness.range.size, data=final2))

plot(kurtosis.range.overlap~kurtosis.range.size, data=final1, xlab="Range size", main="Kurtosis", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(kurtosis.range.overlap~kurtosis.range.size, data=final1))

plot(kurtosis.range.overlap~kurtosis.range.size, data=final2, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
abline(lm(kurtosis.range.overlap~kurtosis.range.size, data=final2))


