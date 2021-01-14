rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

gc()
memory.limit(size=10000000)

library(sna)
library(viridis)

source("codes/functions.R")

dat<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")

dat_NAN<-dat[dat$NAN=="Y",c(10,12)]
rownames(dat_NAN)<-paste(dat[dat$NAN=="Y",3], dat[dat$NAN=="Y",4], sep="_")
dat_NAN[,c(1:2)]<-apply(dat_NAN[,c(1:2)], 2, function(x) as.numeric((x)))

dat_NAN<-na.omit(dat_NAN)

range.sizes.NAN<-numeric()
for(i in 1:nrow(dat_NAN)){
  range.sizes.NAN[i]<-dat_NAN$MAX[i]-dat_NAN$MIN[i]
}
names(range.sizes.NAN)<-rownames(dat_NAN)
range.sizes.NAN<-range.sizes.NAN[range.sizes.NAN>0]

dat_NAN<-dat_NAN[rownames(dat_NAN) %in% names(range.sizes.NAN),]
S_NAN<-nrow(dat_NAN)

overlap.NAN<-overlap.matrix(dat_NAN)

birds<-matrix(nrow=5, ncol=1)
colnames(birds)<-c("NAN")
rownames(birds)<-c("rho", "log.mean.range.size",
                   "log.variance.range.size", 
                   "log.mean.range.overlap",
                   "log.variance.range.overlap")


birds[1,1]<-mean(log(range.sizes.NAN))/(max(dat_NAN$MAX)-min(dat_NAN$MIN))
birds[2,1]<-mean(log(diag(overlap.NAN)))
birds[3,1]<-var(log(diag(overlap.NAN)))

diag(overlap.NAN)<-NA

over.NAN<-overlap.NAN[overlap.NAN > 0]

birds[4,1]<-mean(log(over.NAN), na.rm = TRUE)
birds[5,1]<-var(as.vector(log(over.NAN)), na.rm = TRUE)

meanlog.NAN<-mean(log(range.sizes.NAN))

sdlog.NAN<-sd(log(range.sizes.NAN))

range.sizes <- rlnorm(S_NAN, meanlog=meanlog.NAN, sdlog=sdlog.NAN)

sim.NAN<-list()
for(i in 1:1000){
  sim.NAN[[i]]<-simRangeOver2(meanlog.NAN, sdlog.NAN, S_NAN, truncate=TRUE, 
                              D.min=min(dat_NAN$MIN), D.max=max(dat_NAN$MAX))
}

res_sim<-matrix(nrow=5, ncol=length(sim.NAN))
rownames(res_sim)<-c("rho", "log.mean.range.size", 
                     "log.variance.range.size", 
                     "log.mean.range.overlap", 
                     "log.variance.range.overlap")
sim.NAN.over<-sim.range.size.NAN<-sim.overlap.NAN<-list()
for(i in 1:length(sim.NAN)){
  res_sim[1,i]<-sim.NAN[[i]]$rho
  sim.range.size.NAN[[i]]<-diag(sim.NAN[[i]]$over)
  res_sim[2,i]<-mean(log(diag(sim.NAN[[i]]$over)))
  res_sim[3,i]<-var(log(diag(sim.NAN[[i]]$over)))

  diag(sim.NAN[[i]]$over)<-NA
  
  sim.NAN.over[[i]]<-sim.NAN[[i]]$over[sim.NAN[[i]]$over > 0]
  sim.overlap.NAN[[i]]<-sim.NAN.over[[i]]
  res_sim[4,i]<-mean(log(sim.NAN.over[[i]]), na.rm = TRUE)
  res_sim[5,i]<-var(as.vector(log(sim.NAN.over[[i]])), na.rm = TRUE)
}

hist(sim.range.size.NAN[[1]], col=rgb(0.1, 0.1, 0.1, 0.1), breaks=40, border=F, ylim=c(0,300))
for(i in 2:1000){
  hist(sim.range.size.NAN[[i]], border=F, breaks=40, col=rgb(0.1, 0.1, 0.1, 0.1), add=T)
}

hist(range.sizes.NAN, col=magma(20)[14], breaks=40, border=F, add=T)

median.sim<-numeric()
for(i in 1:1000){
  median.sim[[i]]<-median(sim.overlap.NAN[[i]])
}

hist(median.sim, border=F, col=rgb(0.2, 0.2, 0.2, 0.2))
xx<-median(log(over.NAN))
abline(v=xx, col=magma(20)[14], lwd=2)

hist(sim.overlap.NAN[[1]], col=rgb(0.1, 0.1, 0.1, 0.1), breaks=40, border=F, ylim=c(0,60000))
for(i in 2:1000){
  hist(sim.overlap.NAN[[i]], border=F, breaks=40, col=rgb(0.1, 0.1, 0.1, 0.1), add=T)
}

over.NAN<-na.omit(over.NAN)
hist(log(over.NAN), col=magma(20)[14], breaks=40, border=F, add=T)


layout(matrix(1:5, ncol=5, byrow = T))
par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

hist_rho<-hist(res_sim[1,], breaks=seq(from=min(res_sim[1,])-0.2, 
                                     to=max(res_sim[1,])+0.2, by=0.02), plot=F)
col<-ifelse(hist_rho$breaks < quantile(res_sim[1,], probs=0.025), 
                 rgb(0.1,0.1,0.1,0.1) , ifelse (hist_rho$breaks >= quantile(res_sim[1,], probs=0.975), 
                                                rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist_rho, col=col, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,1], col=magma(20)[14], lwd=2)

hist2<-hist(res_sim[2,], breaks=seq(from=min(res_sim[2,])-0.1, 
                                       to=max(res_sim[2,])+0.1, by=0.01), plot=F)
col<-ifelse(hist2$breaks < quantile(res_sim[2,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist2$breaks >= quantile(res_sim[2,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2, col=col, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,1], col=magma(20)[14], lwd=2)

hist3<-hist(res_sim[3,], breaks=seq(from=min(res_sim[3,])-0.1, 
                                    to=max(res_sim[3,])+0.1, by=0.01), plot=F)
col<-ifelse(hist3$breaks < quantile(res_sim[3,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist3$breaks >= quantile(res_sim[3,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3, col=col, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,1], col=magma(20)[14], lwd=2)

hist4<-hist(res_sim[4,], breaks=seq(from=min(res_sim[4,])-0.1, 
                                    to=max(res_sim[4,])+0.2, by=0.01), plot=F)
col<-ifelse(hist4$breaks < quantile(res_sim[4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist4$breaks >= quantile(res_sim[4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,1], col=magma(20)[14], lwd=2)


hist(res_sim[[1]][5,], col="white", border="gray", xlab="", main="log(variance RO)", breaks=seq(from=min(res_sim[[1]][5,])-1, to=max(res_sim[[1]][5,])+0.1, by=0.02))
sim1_varro_ci<-res_sim[[1]][5,][res_sim[[1]][5,] >= quantile(res_sim[[1]][5,], probs=0.025)]
sim1_varro_ci<-sim1_varro_ci[sim1_varro_ci <= quantile(res_sim[[1]][5,], probs=0.975)]
hist(sim1_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][5,])-0.1, to=max(res_sim[[1]][5,])+1, by=0.02))
abline(v=birds[5,1], col="red")


dev.off()

layout(matrix(1:2, ncol=2))
rho<-c(res_sim[1,], res_sim[[2]][1,], res_sim[[3]][1,], res_sim[[4]][1,], 
       res_sim[[5]][1,])
meanRO<-c(res_sim[[1]][4,], res_sim[[2]][4,], res_sim[[3]][4,], res_sim[[4]][4,], 
          res_sim[[5]][4,])
meanRS<-c(res_sim[[1]][2,], res_sim[[2]][2,], res_sim[[3]][2,], res_sim[[4]][2,], 
          res_sim[[5]][2,])
meanRO<-c(res_sim[[1]][4,], res_sim[[2]][4,], res_sim[[3]][4,], res_sim[[4]][4,], 
          res_sim[[5]][4,])
plot(res_sim[1,], res_sim[4,], ylim=c(6,8), xlim=c(0,0.1))
points(birds[1,1], birds[4,1], pch=16, col="red")
points(birds[2,2], birds[4,2], pch=16, col="red")
points(birds[2,3], birds[4,3], pch=16, col="red")
points(birds[2,4], birds[4,4], pch=16, col="red")
points(birds[2,5], birds[4,5], pch=16, col="red")



net.birds<-list()
net.birds$AMN<-matrix(nrow=S_AMN, ncol=3)
net.birds$AMS<-matrix(nrow=S_AMS, ncol=3)
net.birds$CAN<-matrix(nrow=S_CAN, ncol=3)
net.birds$NAN<-matrix(nrow=S_NAN, ncol=3)
net.birds$GCS<-matrix(nrow=S_GCS, ncol=3)

overlap.AMN[overlap.AMN > 0]<-overlap.AMS[overlap.AMS > 0]<-
  overlap.CAN[overlap.CAN > 0]<-overlap.NAN[overlap.NAN > 0]<-
  overlap.GCS[overlap.GCS > 0]<-1

net.birds$AMN[,1]<-sort(degree(overlap.AMN, gmode="graph"), decreasing=T)
net.birds$AMN[,2]<-sort(betweenness(overlap.AMN, gmode="graph"), decreasing=T)
net.birds$AMN[,3]<-sort(closeness(overlap.AMN, gmode="graph"), decreasing=T)

net.birds$AMS[,1]<-sort(degree(overlap.AMS, gmode="graph"), decreasing=T)
net.birds$AMS[,2]<-sort(betweenness(overlap.AMS, gmode="graph"), decreasing=T)
net.birds$AMS[,3]<-sort(closeness(overlap.AMS, gmode="graph"), decreasing=T)

net.birds$CAN[,1]<-sort(degree(overlap.CAN, gmode="graph"), decreasing=T)
net.birds$CAN[,2]<-sort(betweenness(overlap.CAN, gmode="graph"), decreasing=T)
net.birds$CAN[,3]<-sort(closeness(overlap.CAN, gmode="graph"), decreasing=T)

net.birds$NAN[,1]<-sort(degree(overlap.NAN, gmode="graph"), decreasing=T)
net.birds$NAN[,2]<-sort(betweenness(overlap.NAN, gmode="graph"), decreasing=T)
net.birds$NAN[,3]<-sort(closeness(overlap.NAN, gmode="graph"), decreasing=T)

net.birds$GCS[,1]<-sort(degree(overlap.GCS, gmode="graph"), decreasing=T)
net.birds$GCS[,2]<-sort(betweenness(overlap.GCS, gmode="graph"), decreasing=T)
net.birds$GCS[,3]<-sort(closeness(overlap.GCS, gmode="graph"), decreasing=T)

net.sim<-list()
for(i in 1:1000){
  net.sim[[i]]<-list()
  
  net.sim[[i]]$AMN<-matrix(nrow=S_AMN, ncol=3)
  net.sim[[i]]$AMS<-matrix(nrow=S_AMS, ncol=3)
  net.sim[[i]]$CAN<-matrix(nrow=S_CAN, ncol=3)
  net.sim[[i]]$NAN<-matrix(nrow=S_NAN, ncol=3)
  net.sim[[i]]$GCS<-matrix(nrow=S_GCS, ncol=3)
  
  sim.AMN[[i]]$over[sim.AMN[[i]]$over > 0]<-sim.AMS[[i]]$over[sim.AMS[[i]]$over > 0]<- 
    sim.CAN[[i]]$over[sim.CAN[[i]]$over > 0]<-sim.NAN[[i]]$over[sim.NAN[[i]]$over > 0]<-
    sim.GCS[[i]]$over[sim.GCS[[i]]$over > 0]<-1
  
  net.sim[[i]]$AMN[,1]<-sort(degree(sim.AMN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$AMN[,2]<-sort(betweenness(sim.AMN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$AMN[,3]<-sort(closeness(sim.AMN[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$AMS[,1]<-sort(degree(sim.AMS[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$AMS[,2]<-sort(betweenness(sim.AMS[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$AMS[,3]<-sort(closeness(sim.AMS[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$CAN[,1]<-sort(degree(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CAN[,2]<-sort(betweenness(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CAN[,3]<-sort(closeness(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$NAN[,1]<-sort(degree(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$NAN[,2]<-sort(betweenness(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$NAN[,3]<-sort(closeness(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$GCS[,1]<-sort(degree(sim.GCS[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$GCS[,2]<-sort(betweenness(sim.GCS[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$GCS[,3]<-sort(closeness(sim.GCS[[i]]$over, gmode="graph"), decreasing=T)
}

pdf("figures/oi.pdf")
layout(matrix(1:15, ncol=3, byrow=TRUE))
par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col <- rbPal1(5)[as.numeric(cut(1:5, breaks = 5))]

plot(net.sim[[1]]$AMN[,1], main="Degree", ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMN[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMN[,1], col="gray")
}
lines(net.birds$AMN[,1], col=col[1]) 

plot(net.sim[[1]]$AMN[,2], main="Betweenness", ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMN[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMN[,2], col="gray")
}
lines(net.birds$AMN[,2], col=col[1]) 

plot(net.sim[[1]]$AMN[,3], main="Closeness", ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMN[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMN[,3], col="gray")
}
lines(net.birds$AMN[,3], col=col[1]) 


plot(net.sim[[1]]$AMS[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMS[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMS[,1], col="gray")
}
lines(net.birds$AMS[,1], col=col[2]) 

plot(net.sim[[1]]$AMS[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMS[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMS[,2], col="gray")
}
lines(net.birds$AMS[,2], col=col[2]) 

plot(net.sim[[1]]$AMS[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$AMS[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$AMS[,3], col="gray")
}
lines(net.birds$AMS[,3], col=col[2]) 


plot(net.sim[[1]]$CAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$CAN[,1], col="gray")
}
lines(net.birds$CAN[,1], col=col[3]) 

plot(net.sim[[1]]$CAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$CAN[,2], col="gray")
}
lines(net.birds$CAN[,2], col=col[3]) 

plot(net.sim[[1]]$CAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$CAN[,3], col="gray")
}
lines(net.birds$CAN[,3], col=col[3]) 


plot(net.sim[[1]]$NAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$NAN[,1], col="gray")
}
lines(net.birds$NAN[,1], col=col[4]) 

plot(net.sim[[1]]$NAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$NAN[,2], col="gray")
}
lines(net.birds$NAN[,2], col=col[4]) 

plot(net.sim[[1]]$NAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$NAN[,3], col="gray")
}
lines(net.birds$NAN[,3], col=col[4]) 


plot(net.sim[[1]]$GCS[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$GCS[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$GCS[,1], col="gray")
}
lines(net.birds$GCS[,1], col=col[4]) 

plot(net.sim[[1]]$GCS[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$GCS[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$GCS[,2], col="gray")
}
lines(net.birds$GCS[,2], col=col[4]) 

plot(net.sim[[1]]$GCS[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$GCS[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$GCS[,3], col="gray")
}
lines(net.birds$GCS[,3], col=col[4]) 

dev.off()