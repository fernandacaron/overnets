rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

meanlog<-seq(from = -5, to = 1, by = 0.012)
S<-1000

final1<-final2<-as.data.frame(matrix(ncol = 10, nrow = length(meanlog)))
colnames(final1)<-colnames(final2)<-c("rho", "mean.range.size", "log.mean.range.size", 
                                      "mean.range.overlap", "log.mean.range.overlap",
                                      "variance.range.size", "log.variance.range.size",
                                      "variance.range.overlap", "log.variance.range.overlap",
                                      "prop.zeros.overlap")

for (i in 1:length(meanlog)){
  progress(i, max.value = length(meanlog))
  dat<-simRangeOver(meanlog=meanlog[i], S=S, truncate=TRUE)
  overlap<-dat$over
  rho<-dat$rho
  
  #range size stats
  mean.range.size<-mean(diag(overlap))
  log.mean.range.size<-mean(log(diag(overlap)))
  variance.range.size<-var(diag(overlap))
  log.variance.range.size<-var(log(diag(overlap)))
  
  #range overlap stats
  diag(overlap)<-NA
  x<-overlap
  x[x > 0]<-1
  prop.zeros.overlap<-mean(x[upper.tri(x)])
  overlap<-overlap[overlap > 0]
  mean.range.overlap<-mean(overlap, na.rm = TRUE)
  log.mean.range.overlap<-mean(log(overlap), na.rm = TRUE)
  variance.range.overlap<-var(as.vector(overlap), na.rm = TRUE)
  log.variance.range.overlap<-var(as.vector(log(overlap)), na.rm = TRUE)

  
  res<-cbind(rho, mean.range.size, log.mean.range.size, mean.range.overlap, 
             log.mean.range.overlap, variance.range.size, log.variance.range.size, 
             variance.range.overlap, log.variance.range.overlap, prop.zeros.overlap)
  
  final1[i,]<-res
}

for (i in 1:length(meanlog)){
  progress(i, max.value = length(meanlog))
  
  dat<-simRangeOver(meanlog=meanlog[i], S=S, truncate=FALSE)
  overlap<-dat$over
  rho<-dat$rho
  
  #range size stats
  mean.range.size<-mean(diag(overlap))
  log.mean.range.size<-mean(log(diag(overlap)))
  variance.range.size<-var(diag(overlap))
  log.variance.range.size<-var(log(diag(overlap)))

  #range overlap stats
  diag(overlap)<-NA
  x<-overlap
  x[x > 0]<-1
  prop.zeros.overlap<-mean(x[upper.tri(x)])
  overlap <- overlap[overlap > 0]
  mean.range.overlap<-mean(overlap, na.rm = TRUE)
  log.mean.range.overlap<-mean(log(overlap), na.rm = TRUE)
  variance.range.overlap<-var(as.vector(overlap), na.rm = TRUE)
  log.variance.range.overlap<-var(as.vector(log(overlap)), na.rm = TRUE)

  
  res<-cbind(rho, mean.range.size, log.mean.range.size, mean.range.overlap, 
             log.mean.range.overlap, variance.range.size, log.variance.range.size, 
             variance.range.overlap, log.variance.range.overlap, prop.zeros.overlap)
  
  final2[i,]<-res
}

#Proportion 0s
pdf("figures/Figure2.pdf")
layout(matrix(1:2, ncol=2))

plot(prop.zeros.overlap~rho, data=final1, xlab=expression(rho), main="", ylab = "Proportion of overlapping species", col = rgb(50, 100, 142, 120, maxColorValue = 255), pch = 16)
abline(lm(prop.zeros.overlap~rho, data=final1))

plot(prop.zeros.overlap~rho, data=final2, xlab=expression(rho), main="", ylab = "Proportion of overlapping species", col = rgb(182, 48, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(prop.zeros.overlap~rho, data=final2))

dev.off()

#Mean and Variance
pdf("figures/Figure3.pdf")
layout(matrix(1:4, ncol=2))

plot(log.mean.range.overlap~log.mean.range.size, data=final1, xlab=expression(paste("log(", rho, ")")), main="Mean", ylab = "log(Range overlap)", col = rgb(50, 100, 142, 120, maxColorValue = 255), pch = 16)
abline(lm(log.mean.range.overlap~log.mean.range.size, data=final1))

plot(log.mean.range.overlap~log.mean.range.size, data=final2, xlab=expression(paste("log(", rho, ")")), main="", ylab = "log(Range overlap)", col = rgb(182, 48, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(log.mean.range.overlap~log.mean.range.size, data=final2))

plot(log.variance.range.overlap~log.variance.range.size, data=final1, xlab=expression(paste("log(", rho, ")")), main="Variance", ylab = "log(Range overlap)", col = rgb(50, 100, 142, 120, maxColorValue = 255), pch = 16)
abline(lm(log.variance.range.overlap~log.variance.range.size, data=final1))

plot(log.variance.range.overlap~log.variance.range.size, data=final2, xlab=expression(paste("log(", rho, ")")), main="", ylab = "log(Range overlap)", col = rgb(182, 48, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(log.variance.range.overlap~log.variance.range.size, data=final2))

dev.off()

pdf("figures/FigureS1.pdf")
layout(matrix(1:4, ncol=2))

plot(mean.range.overlap~mean.range.size, data=final1, xlab=expression(rho), main="Mean", ylab = "Range overlap", col = rgb(50, 100, 142, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final1))

plot(mean.range.overlap~mean.range.size, data=final2, xlab=expression(rho), main="", ylab = "Range overlap", col = rgb(182, 48, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.range.overlap~mean.range.size, data=final2))

plot(variance.range.overlap~variance.range.size, data=final1, xlab=expression(rho), main="Variance", ylab = "Range overlap", col = rgb(50, 100, 142, 120, maxColorValue = 255), pch = 16)
abline(lm(variance.range.overlap~variance.range.size, data=final1))

plot(variance.range.overlap~variance.range.size, data=final2, xlab=expression(rho), main="", ylab = "Range overlap", col = rgb(182, 48, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(variance.range.overlap~variance.range.size, data=final2))

dev.off()
