rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

gc()
memory.limit(size=10000000)

library(sna)
library(viridis)

source("codes/functions.R")

dat<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")

dat_MAH<-dat[dat$MAH=="Y",c(10,12)]
rownames(dat_MAH)<-paste(dat[dat$MAH=="Y",3], dat[dat$MAH=="Y",4], sep="_")
dat_MAH[,c(1:2)]<-apply(dat_MAH[,c(1:2)], 2, function(x) as.numeric((x)))

dat_CDH<-dat[dat$CDH=="Y",c(10,12)]
rownames(dat_CDH)<-paste(dat[dat$CDH=="Y",3], dat[dat$CDH=="Y",4], sep="_")
dat_CDH[,c(1:2)]<-apply(dat_CDH[,c(1:2)], 2, function(x) as.numeric((x)))

dat_CAN<-dat[dat$CAN=="Y",c(10,12)]
rownames(dat_CAN)<-paste(dat[dat$CAN=="Y",3], dat[dat$CAN=="Y",4], sep="_")
dat_CAN[,c(1:2)]<-apply(dat_CAN[,c(1:2)], 2, function(x) as.numeric((x)))

dat_NAN<-dat[dat$NAN=="Y",c(10,12)]
rownames(dat_NAN)<-paste(dat[dat$NAN=="Y",3], dat[dat$NAN=="Y",4], sep="_")
dat_NAN[,c(1:2)]<-apply(dat_NAN[,c(1:2)], 2, function(x) as.numeric((x)))

dat_TEP<-dat[dat$TEP=="Y",c(10,12)]
rownames(dat_TEP)<-paste(dat[dat$TEP=="Y",3], dat[dat$TEP=="Y",4], sep="_")
dat_TEP[,c(1:2)]<-apply(dat_TEP[,c(1:2)], 2, function(x) as.numeric((x)))

dat_MAH<-na.omit(dat_MAH)
dat_CDH<-na.omit(dat_CDH)
dat_CAN<-na.omit(dat_CAN)
dat_NAN<-na.omit(dat_NAN)
dat_TEP<-na.omit(dat_TEP)

range.sizes.MAH<-numeric()
for(i in 1:nrow(dat_MAH)){
  range.sizes.MAH[i]<-dat_MAH$MAX[i]-dat_MAH$MIN[i]
}
names(range.sizes.MAH)<-rownames(dat_MAH)
range.sizes.MAH<-range.sizes.MAH[range.sizes.MAH>0]

range.sizes.CDH<-numeric()
for(i in 1:nrow(dat_CDH)){
  range.sizes.CDH[i]<-dat_CDH$MAX[i]-dat_CDH$MIN[i]
}
names(range.sizes.CDH)<-rownames(dat_CDH)
range.sizes.CDH<-range.sizes.CDH[range.sizes.CDH>0]

range.sizes.CAN<-numeric()
for(i in 1:nrow(dat_CAN)){
  range.sizes.CAN[i]<-dat_CAN$MAX[i]-dat_CAN$MIN[i]
}
names(range.sizes.CAN)<-rownames(dat_CAN)
range.sizes.CAN<-range.sizes.CAN[range.sizes.CAN>0]

range.sizes.NAN<-numeric()
for(i in 1:nrow(dat_NAN)){
  range.sizes.NAN[i]<-dat_NAN$MAX[i]-dat_NAN$MIN[i]
}
names(range.sizes.NAN)<-rownames(dat_NAN)
range.sizes.NAN<-range.sizes.NAN[range.sizes.NAN>0]

range.sizes.TEP<-numeric()
for(i in 1:nrow(dat_TEP)){
  range.sizes.TEP[i]<-dat_TEP$MAX[i]-dat_TEP$MIN[i]
}
names(range.sizes.TEP)<-rownames(dat_TEP)
range.sizes.TEP<-range.sizes.TEP[range.sizes.TEP>0]

dat_MAH<-dat_MAH[rownames(dat_MAH) %in% names(range.sizes.MAH),]
dat_CDH<-dat_CDH[rownames(dat_CDH) %in% names(range.sizes.CDH),]
dat_CAN<-dat_CAN[rownames(dat_CAN) %in% names(range.sizes.CAN),]
dat_NAN<-dat_NAN[rownames(dat_NAN) %in% names(range.sizes.NAN),]
dat_TEP<-dat_TEP[rownames(dat_TEP) %in% names(range.sizes.TEP),]

S_MAH<-nrow(dat_MAH)
S_CDH<-nrow(dat_CDH)
S_CAN<-nrow(dat_CAN)
S_NAN<-nrow(dat_NAN)
S_TEP<-nrow(dat_TEP)

overlap.MAH<-overlap.matrix(dat_MAH)
overlap.CDH<-overlap.matrix(dat_CDH)
overlap.CAN<-overlap.matrix(dat_CAN)
overlap.NAN<-overlap.matrix(dat_NAN)
overlap.TEP<-overlap.matrix(dat_TEP)

birds<-matrix(nrow=5, ncol=5)
colnames(birds)<-c("MAH", "CDH", "CAN", "NAN", "TEP")
rownames(birds)<-c("rho", "log.mean.range.size",
                   "log.variance.range.size", 
                   "log.mean.range.overlap",
                   "log.variance.range.overlap")


birds[1,1]<-mean(log(range.sizes.MAH))/(max(dat_MAH$MAX)-min(dat_MAH$MIN))
birds[1,2]<-mean(log(range.sizes.CDH))/(max(dat_CDH$MAX)-min(dat_CDH$MIN))
birds[1,3]<-mean(log(range.sizes.CAN))/(max(dat_CAN$MAX)-min(dat_CAN$MIN))
birds[1,4]<-mean(log(range.sizes.NAN))/(max(dat_NAN$MAX)-min(dat_NAN$MIN))
birds[1,5]<-mean(log(range.sizes.TEP))/(max(dat_TEP$MAX)-min(dat_TEP$MIN))

birds[2,1]<-mean(log(diag(overlap.MAH)))
birds[2,2]<-mean(log(diag(overlap.CDH)))
birds[2,3]<-mean(log(diag(overlap.CAN)))
birds[2,4]<-mean(log(diag(overlap.NAN)))
birds[2,5]<-mean(log(diag(overlap.TEP)))

birds[3,1]<-var(log(diag(overlap.MAH)))
birds[3,2]<-var(log(diag(overlap.CDH)))
birds[3,3]<-var(log(diag(overlap.CAN)))
birds[3,4]<-var(log(diag(overlap.NAN)))
birds[3,5]<-var(log(diag(overlap.TEP)))

over.MAH<-overlap.MAH
over.CDH<-overlap.CDH
over.CAN<-overlap.CAN
over.NAN<-overlap.NAN
over.TEP<-overlap.TEP

diag(over.MAH)<-diag(over.CDH)<-diag(over.CAN)<-diag(over.NAN)<-diag(over.TEP)<-NA

over.MAH<-over.MAH[over.MAH > 0]
over.CDH<-over.CDH[over.CDH > 0]
over.CAN<-over.CAN[over.CAN > 0]
over.NAN<-over.NAN[over.NAN > 0]
over.TEP<-over.TEP[over.TEP > 0]

birds[4,1]<-mean(log(over.MAH), na.rm = TRUE)
birds[4,2]<-mean(log(over.CDH), na.rm = TRUE)
birds[4,3]<-mean(log(over.CAN), na.rm = TRUE)
birds[4,4]<-mean(log(over.NAN), na.rm = TRUE)
birds[4,5]<-mean(log(over.TEP), na.rm = TRUE)

birds[5,1]<-var(as.vector(log(over.MAH)), na.rm = TRUE)
birds[5,2]<-var(as.vector(log(over.CDH)), na.rm = TRUE)
birds[5,3]<-var(as.vector(log(over.CAN)), na.rm = TRUE)
birds[5,4]<-var(as.vector(log(over.NAN)), na.rm = TRUE)
birds[5,5]<-var(as.vector(log(over.TEP)), na.rm = TRUE)

meanlog.MAH<-mean(log(range.sizes.MAH))
meanlog.CDH<-mean(log(range.sizes.CDH))
meanlog.CAN<-mean(log(range.sizes.CAN))
meanlog.NAN<-mean(log(range.sizes.NAN))
meanlog.TEP<-mean(log(range.sizes.TEP))

sdlog.MAH<-sd(log(range.sizes.MAH))
sdlog.CDH<-sd(log(range.sizes.CDH))
sdlog.CAN<-sd(log(range.sizes.CAN))
sdlog.NAN<-sd(log(range.sizes.NAN))
sdlog.TEP<-sd(log(range.sizes.TEP))

sim.MAH<-sim.CDH<-sim.CAN<-sim.NAN<-sim.TEP<-list()
for(i in 1:1000){
  sim.MAH[[i]]<-simRangeOver2(meanlog.MAH, sdlog.MAH, S_MAH, truncate=TRUE, 
                              D.min=min(dat_MAH$MIN), D.max=max(dat_MAH$MAX))
  sim.CDH[[i]]<-simRangeOver2(meanlog.CDH, sdlog.CDH, S_CDH, truncate=TRUE, 
                              D.min=min(dat_CDH$MIN), D.max=max(dat_CDH$MAX))
  sim.CAN[[i]]<-simRangeOver2(meanlog.CAN, sdlog.CAN, S_CAN, truncate=TRUE, 
                              D.min=min(dat_CAN$MIN), D.max=max(dat_CAN$MAX))
  sim.NAN[[i]]<-simRangeOver2(meanlog.NAN, sdlog.NAN, S_NAN, truncate=TRUE, 
                              D.min=min(dat_NAN$MIN), D.max=max(dat_NAN$MAX))
  sim.TEP[[i]]<-simRangeOver2(meanlog.TEP, sdlog.TEP, S_TEP, truncate=TRUE, 
                              D.min=min(dat_TEP$MIN), D.max=max(dat_TEP$MAX))
}

res_sim<-list()
res_sim[[1]]<-matrix(nrow=5, ncol=length(sim.MAH))
res_sim[[2]]<-matrix(nrow=5, ncol=length(sim.CDH))
res_sim[[3]]<-matrix(nrow=5, ncol=length(sim.CAN))
res_sim[[4]]<-matrix(nrow=5, ncol=length(sim.NAN))
res_sim[[5]]<-matrix(nrow=5, ncol=length(sim.TEP))
rownames(res_sim[[1]])<-rownames(res_sim[[2]])<-rownames(res_sim[[4]])<-
  rownames(res_sim[[5]])<-rownames(res_sim[[3]])<-c("rho", "log.mean.range.size",
                                                    "log.variance.range.size", 
                                                    "log.mean.range.overlap", 
                                                    "log.variance.range.overlap")

sim.MAH.over<-sim.CDH.over<-sim.CAN.over<-sim.NAN.over<-sim.TEP.over<-list()
for(i in 1:length(sim.NAN)){
  res_sim[[1]][1,i]<-sim.MAH[[i]]$rho
  res_sim[[2]][1,i]<-sim.CDH[[i]]$rho
  res_sim[[3]][1,i]<-sim.CAN[[i]]$rho
  res_sim[[4]][1,i]<-sim.NAN[[i]]$rho
  res_sim[[5]][1,i]<-sim.TEP[[i]]$rho

  res_sim[[1]][2,i]<-mean(log(diag(sim.MAH[[i]]$over)))
  res_sim[[2]][2,i]<-mean(log(diag(sim.CDH[[i]]$over)))
  res_sim[[3]][2,i]<-mean(log(diag(sim.CAN[[i]]$over)))
  res_sim[[4]][2,i]<-mean(log(diag(sim.NAN[[i]]$over)))
  res_sim[[5]][2,i]<-mean(log(diag(sim.TEP[[i]]$over)))
  
  res_sim[[1]][3,i]<-var(log(diag(sim.MAH[[i]]$over)))
  res_sim[[2]][3,i]<-var(log(diag(sim.CDH[[i]]$over)))
  res_sim[[3]][3,i]<-var(log(diag(sim.CAN[[i]]$over)))
  res_sim[[4]][3,i]<-var(log(diag(sim.NAN[[i]]$over)))
  res_sim[[5]][3,i]<-var(log(diag(sim.TEP[[i]]$over)))
  
  sim.MAH.over[[i]]<-sim.MAH[[i]]$over
  sim.CDH.over[[i]]<-sim.CDH[[i]]$over
  sim.CAN.over[[i]]<-sim.CAN[[i]]$over
  sim.NAN.over[[i]]<-sim.NAN[[i]]$over
  sim.TEP.over[[i]]<-sim.TEP[[i]]$over
  
  diag(sim.MAH.over[[i]])<-NA
  diag(sim.CDH.over[[i]])<-NA
  diag(sim.CAN.over[[i]])<-NA
  diag(sim.NAN.over[[i]])<-NA
  diag(sim.TEP.over[[i]])<-NA
  
  sim.MAH.over[[i]]<-sim.MAH.over[[i]][sim.MAH.over[[i]] > 0]
  sim.CDH.over[[i]]<-sim.CDH.over[[i]][sim.CDH.over[[i]] > 0]
  sim.CAN.over[[i]]<-sim.CAN.over[[i]][sim.CAN.over[[i]] > 0]
  sim.NAN.over[[i]]<-sim.NAN.over[[i]][sim.NAN.over[[i]] > 0]
  sim.TEP.over[[i]]<-sim.TEP.over[[i]][sim.TEP.over[[i]] > 0]

  res_sim[[1]][4,i]<-mean(log(sim.MAH.over[[i]]), na.rm = TRUE)
  res_sim[[2]][4,i]<-mean(log(sim.CDH.over[[i]]), na.rm = TRUE)
  res_sim[[3]][4,i]<-mean(log(sim.CAN.over[[i]]), na.rm = TRUE)
  res_sim[[4]][4,i]<-mean(log(sim.NAN.over[[i]]), na.rm = TRUE)
  res_sim[[5]][4,i]<-mean(log(sim.TEP.over[[i]]), na.rm = TRUE)
  
  res_sim[[1]][5,i]<-var(as.vector(log(sim.MAH.over[[i]])), na.rm = TRUE)
  res_sim[[2]][5,i]<-var(as.vector(log(sim.CDH.over[[i]])), na.rm = TRUE)
  res_sim[[3]][5,i]<-var(as.vector(log(sim.CAN.over[[i]])), na.rm = TRUE)
  res_sim[[4]][5,i]<-var(as.vector(log(sim.NAN.over[[i]])), na.rm = TRUE)
  res_sim[[5]][5,i]<-var(as.vector(log(sim.TEP.over[[i]])), na.rm = TRUE)
}

load(file="data/res_sim_NT.RData")
pdf("figures/Figure6A.pdf")
layout(matrix(1:25, ncol=5, byrow = T))
par(mar=c(4, 4, 0.5, 0.5))

hist1_1<-hist(res_sim[[1]][1,], breaks=seq(from=min(res_sim[[1]][1,])-0.00001, 
                                       to=max(res_sim[[1]][1,])+0.00001, by=0.0000015), plot=F)
col1_1<-ifelse(hist1_1$breaks < quantile(res_sim[[1]][1,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1), ifelse (hist1_1$breaks >= 
                                            quantile(res_sim[[1]][1,], probs=0.975), 
                                          rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist1_1, col=col1_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,1], col=magma(20)[14], lwd=2)

hist1_2<-hist(res_sim[[1]][2,], breaks=seq(from=min(res_sim[[1]][2,])-0.1,
                                           to=max(res_sim[[1]][2,])+0.1, by=0.01), plot=F)
col1_2<-ifelse(hist1_2$breaks < quantile(res_sim[[1]][2,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_2$breaks >= quantile(res_sim[[1]][2,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_2, col=col1_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,1], col=magma(20)[14], lwd=2)

hist1_3<-hist(res_sim[[1]][3,], breaks=seq(from=min(res_sim[[1]][3,])-0.1, 
                                           to=max(res_sim[[1]][3,])+0.1, by=0.008), plot=F)
col1_3<-ifelse(hist1_3$breaks < quantile(res_sim[[1]][3,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_3$breaks >= quantile(res_sim[[1]][3,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_3, col=col1_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,1], col=magma(20)[14], lwd=2)

hist1_4<-hist(res_sim[[1]][4,], breaks=seq(from=min(res_sim[[1]][4,])-0.1, 
                                           to=max(res_sim[[1]][4,])+0.3, by=0.02), plot=F)
col<-ifelse(hist1_4$breaks < quantile(res_sim[[1]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_4$breaks >= quantile(res_sim[[1]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,1], col=magma(20)[14], lwd=2)

hist1_5<-hist(res_sim[[1]][5,], breaks=seq(from=min(res_sim[[1]][5,])-0.5, 
                                           to=max(res_sim[[1]][5,])+0.1, by=0.02), plot=F)
col1_5<-ifelse(hist1_5$breaks < quantile(res_sim[[1]][5,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_5$breaks >= quantile(res_sim[[1]][5,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_5, col=col1_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,1], col=magma(20)[14], lwd=2)



hist2_1<-hist(res_sim[[2]][1,], breaks=seq(from=min(res_sim[[2]][1,])-0.00001, 
                                           to=max(res_sim[[2]][1,])+0.00001, by=0.000002), plot=F)
col2_1<-ifelse(hist2_1$breaks < quantile(res_sim[[2]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist2_1$breaks >= 
                                               quantile(res_sim[[2]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist2_1, col=col2_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,2], col=magma(20)[14], lwd=2)

hist2_2<-hist(res_sim[[2]][2,], breaks=seq(from=min(res_sim[[2]][2,])-0.2,
                                           to=max(res_sim[[2]][2,])+0.1, by=0.01), plot=F)
col2_2<-ifelse(hist2_2$breaks < quantile(res_sim[[2]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_2$breaks >= quantile(res_sim[[2]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_2, col=col2_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,2], col=magma(20)[14], lwd=2)

hist2_3<-hist(res_sim[[2]][3,], breaks=seq(from=min(res_sim[[2]][3,])-0.1, 
                                           to=max(res_sim[[2]][3,])+0.1, by=0.008), plot=F)
col2_3<-ifelse(hist2_3$breaks < quantile(res_sim[[2]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_3$breaks >= quantile(res_sim[[2]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_3, col=col2_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,2], col=magma(20)[14], lwd=2)

hist2_4<-hist(res_sim[[2]][4,], breaks=seq(from=min(res_sim[[2]][4,])-0.1, 
                                           to=max(res_sim[[2]][4,])+0.3, by=0.02), plot=F)
col<-ifelse(hist2_4$breaks < quantile(res_sim[[2]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_4$breaks >= quantile(res_sim[[2]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,2], col=magma(20)[14], lwd=2)

hist2_5<-hist(res_sim[[2]][5,], breaks=seq(from=min(res_sim[[2]][5,])-0.5, 
                                           to=max(res_sim[[2]][5,])+0.1, by=0.02), plot=F)
col2_5<-ifelse(hist2_5$breaks < quantile(res_sim[[2]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_5$breaks >= quantile(res_sim[[2]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_5, col=col2_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,2], col=magma(20)[14], lwd=2)




hist3_1<-hist(res_sim[[3]][1,], breaks=seq(from=min(res_sim[[3]][1,])-0.00001, 
                                           to=max(res_sim[[3]][1,])+0.00001, by=0.0000015), plot=F)
col3_1<-ifelse(hist3_1$breaks < quantile(res_sim[[3]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist3_1$breaks >= 
                                               quantile(res_sim[[3]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist3_1, col=col3_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,3], col=magma(20)[14], lwd=2)

hist3_2<-hist(res_sim[[3]][2,], breaks=seq(from=min(res_sim[[3]][2,])-0.1,
                                           to=max(res_sim[[3]][2,])+0.1, by=0.01), plot=F)
col3_2<-ifelse(hist3_2$breaks < quantile(res_sim[[3]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_2$breaks >= quantile(res_sim[[3]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_2, col=col3_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,3], col=magma(20)[14], lwd=2)

hist3_3<-hist(res_sim[[3]][3,], breaks=seq(from=min(res_sim[[3]][3,])-0.1, 
                                           to=max(res_sim[[3]][3,])+0.1, by=0.008), plot=F)
col3_3<-ifelse(hist3_3$breaks < quantile(res_sim[[3]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_3$breaks >= quantile(res_sim[[3]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_3, col=col3_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,3], col=magma(20)[14], lwd=2)

hist3_4<-hist(res_sim[[3]][4,], breaks=seq(from=min(res_sim[[3]][4,])-0.1, 
                                           to=max(res_sim[[3]][4,])+0.2, by=0.01), plot=F)
col<-ifelse(hist3_4$breaks < quantile(res_sim[[3]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_4$breaks >= quantile(res_sim[[3]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,3], col=magma(20)[14], lwd=2)

hist3_5<-hist(res_sim[[3]][5,], breaks=seq(from=min(res_sim[[3]][5,])-0.3, 
                                           to=max(res_sim[[3]][5,])+0.1, by=0.01), plot=F)
col3_5<-ifelse(hist3_5$breaks < quantile(res_sim[[3]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_5$breaks >= quantile(res_sim[[3]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_5, col=col3_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,3], col=magma(20)[14], lwd=2)



hist4_1<-hist(res_sim[[4]][1,], breaks=seq(from=min(res_sim[[4]][1,])-0.00001, 
                                           to=max(res_sim[[4]][1,])+0.00001, by=0.0000015), plot=F)
col4_1<-ifelse(hist4_1$breaks < quantile(res_sim[[4]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist4_1$breaks >= 
                                               quantile(res_sim[[4]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist4_1, col=col4_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,4], col=magma(20)[14], lwd=2)

hist4_2<-hist(res_sim[[4]][2,], breaks=seq(from=min(res_sim[[4]][2,])-0.1,
                                           to=max(res_sim[[4]][2,])+0.1, by=0.01), plot=F)
col4_2<-ifelse(hist4_2$breaks < quantile(res_sim[[4]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_2$breaks >= quantile(res_sim[[4]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_2, col=col4_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,4], col=magma(20)[14], lwd=2)

hist4_3<-hist(res_sim[[4]][3,], breaks=seq(from=min(res_sim[[4]][3,])-0.1, 
                                           to=max(res_sim[[4]][3,])+0.1, by=0.008), plot=F)
col4_3<-ifelse(hist4_3$breaks < quantile(res_sim[[4]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_3$breaks >= quantile(res_sim[[4]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_3, col=col4_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,4], col=magma(20)[14], lwd=2)

hist4_4<-hist(res_sim[[4]][4,], breaks=seq(from=min(res_sim[[4]][4,])-0.1, 
                                           to=max(res_sim[[4]][4,])+0.2, by=0.01), plot=F)
col<-ifelse(hist4_4$breaks < quantile(res_sim[[4]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_4$breaks >= quantile(res_sim[[4]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,4], col=magma(20)[14], lwd=2)

hist4_5<-hist(res_sim[[4]][5,], breaks=seq(from=min(res_sim[[4]][5,])-0.4, 
                                           to=max(res_sim[[4]][5,])+0.1, by=0.01), plot=F)
col4_5<-ifelse(hist4_5$breaks < quantile(res_sim[[4]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_5$breaks >= quantile(res_sim[[4]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_5, col=col4_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,4], col=magma(20)[14], lwd=2)



hist5_1<-hist(res_sim[[5]][1,], breaks=seq(from=min(res_sim[[5]][1,])-0.00001, 
                                           to=max(res_sim[[5]][1,])+0.00001, by=0.0000025), plot=F)
col5_1<-ifelse(hist5_1$breaks < quantile(res_sim[[5]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist5_1$breaks >= 
                                               quantile(res_sim[[5]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist5_1, col=col5_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,5], col=magma(20)[14], lwd=2)

hist5_2<-hist(res_sim[[5]][2,], breaks=seq(from=min(res_sim[[5]][2,])-0.1,
                                           to=max(res_sim[[5]][2,])+0.1, by=0.01), plot=F)
col5_2<-ifelse(hist5_2$breaks < quantile(res_sim[[5]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_2$breaks >= quantile(res_sim[[5]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_2, col=col5_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,5], col=magma(20)[14], lwd=2)

hist5_3<-hist(res_sim[[5]][3,], breaks=seq(from=min(res_sim[[5]][3,])-0.1, 
                                           to=max(res_sim[[5]][3,])+0.1, by=0.008), plot=F)
col5_3<-ifelse(hist5_3$breaks < quantile(res_sim[[5]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_3$breaks >= quantile(res_sim[[5]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_3, col=col5_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,5], col=magma(20)[14], lwd=2)

hist5_4<-hist(res_sim[[5]][4,], breaks=seq(from=min(res_sim[[5]][4,])-0.1, 
                                           to=max(res_sim[[5]][4,])+0.2, by=0.02), plot=F)
col<-ifelse(hist5_4$breaks < quantile(res_sim[[5]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_4$breaks >= quantile(res_sim[[5]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,5], col=magma(20)[14], lwd=2)

hist5_5<-hist(res_sim[[5]][5,], breaks=seq(from=min(res_sim[[5]][5,])-0.4, 
                                           to=max(res_sim[[5]][5,])+0.1, by=0.02), plot=F)
col5_5<-ifelse(hist5_5$breaks < quantile(res_sim[[5]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_5$breaks >= quantile(res_sim[[5]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_5, col=col5_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,5], col=magma(20)[14], lwd=2)

dev.off()

load(file="data/res_sim_T.RData")
pdf("figures/Figure6B.pdf")
layout(matrix(1:25, ncol=5, byrow = T))
par(mar=c(4, 4, 0.5, 0.5))

hist1_1<-hist(res_sim[[1]][1,], breaks=seq(from=min(res_sim[[1]][1,])-0.00001, 
                                           to=max(res_sim[[1]][1,])+0.00001, by=0.0000015), plot=F)
col1_1<-ifelse(hist1_1$breaks < quantile(res_sim[[1]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist1_1$breaks >= 
                                               quantile(res_sim[[1]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist1_1, col=col1_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,1], col=magma(20)[14], lwd=2)

hist1_2<-hist(res_sim[[1]][2,], breaks=seq(from=min(res_sim[[1]][2,])-0.1,
                                           to=max(res_sim[[1]][2,])+0.1, by=0.01), plot=F)
col1_2<-ifelse(hist1_2$breaks < quantile(res_sim[[1]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_2$breaks >= quantile(res_sim[[1]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_2, col=col1_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,1], col=magma(20)[14], lwd=2)

hist1_3<-hist(res_sim[[1]][3,], breaks=seq(from=min(res_sim[[1]][3,])-0.1, 
                                           to=max(res_sim[[1]][3,])+0.1, by=0.008), plot=F)
col1_3<-ifelse(hist1_3$breaks < quantile(res_sim[[1]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_3$breaks >= quantile(res_sim[[1]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_3, col=col1_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,1], col=magma(20)[14], lwd=2)

hist1_4<-hist(res_sim[[1]][4,], breaks=seq(from=min(res_sim[[1]][4,])-0.1, 
                                           to=max(res_sim[[1]][4,])+0.3, by=0.02), plot=F)
col<-ifelse(hist1_4$breaks < quantile(res_sim[[1]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_4$breaks >= quantile(res_sim[[1]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,1], col=magma(20)[14], lwd=2)

hist1_5<-hist(res_sim[[1]][5,], breaks=seq(from=min(res_sim[[1]][5,])-0.5, 
                                           to=max(res_sim[[1]][5,])+0.1, by=0.02), plot=F)
col1_5<-ifelse(hist1_5$breaks < quantile(res_sim[[1]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist1_5$breaks >= quantile(res_sim[[1]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist1_5, col=col1_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,1], col=magma(20)[14], lwd=2)



hist2_1<-hist(res_sim[[2]][1,], breaks=seq(from=min(res_sim[[2]][1,])-0.00001, 
                                           to=max(res_sim[[2]][1,])+0.00001, by=0.000002), plot=F)
col2_1<-ifelse(hist2_1$breaks < quantile(res_sim[[2]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist2_1$breaks >= 
                                               quantile(res_sim[[2]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist2_1, col=col2_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,2], col=magma(20)[14], lwd=2)

hist2_2<-hist(res_sim[[2]][2,], breaks=seq(from=min(res_sim[[2]][2,])-0.2,
                                           to=max(res_sim[[2]][2,])+0.1, by=0.01), plot=F)
col2_2<-ifelse(hist2_2$breaks < quantile(res_sim[[2]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_2$breaks >= quantile(res_sim[[2]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_2, col=col2_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,2], col=magma(20)[14], lwd=2)

hist2_3<-hist(res_sim[[2]][3,], breaks=seq(from=min(res_sim[[2]][3,])-0.1, 
                                           to=max(res_sim[[2]][3,])+0.1, by=0.008), plot=F)
col2_3<-ifelse(hist2_3$breaks < quantile(res_sim[[2]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_3$breaks >= quantile(res_sim[[2]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_3, col=col2_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,2], col=magma(20)[14], lwd=2)

hist2_4<-hist(res_sim[[2]][4,], breaks=seq(from=min(res_sim[[2]][4,])-0.1, 
                                           to=max(res_sim[[2]][4,])+0.3, by=0.02), plot=F)
col<-ifelse(hist2_4$breaks < quantile(res_sim[[2]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_4$breaks >= quantile(res_sim[[2]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,2], col=magma(20)[14], lwd=2)

hist2_5<-hist(res_sim[[2]][5,], breaks=seq(from=min(res_sim[[2]][5,])-0.5, 
                                           to=max(res_sim[[2]][5,])+0.1, by=0.02), plot=F)
col2_5<-ifelse(hist2_5$breaks < quantile(res_sim[[2]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist2_5$breaks >= quantile(res_sim[[2]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist2_5, col=col2_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,2], col=magma(20)[14], lwd=2)




hist3_1<-hist(res_sim[[3]][1,], breaks=seq(from=min(res_sim[[3]][1,])-0.00001, 
                                           to=max(res_sim[[3]][1,])+0.00001, by=0.0000015), plot=F)
col3_1<-ifelse(hist3_1$breaks < quantile(res_sim[[3]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist3_1$breaks >= 
                                               quantile(res_sim[[3]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist3_1, col=col3_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,3], col=magma(20)[14], lwd=2)

hist3_2<-hist(res_sim[[3]][2,], breaks=seq(from=min(res_sim[[3]][2,])-0.1,
                                           to=max(res_sim[[3]][2,])+0.1, by=0.01), plot=F)
col3_2<-ifelse(hist3_2$breaks < quantile(res_sim[[3]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_2$breaks >= quantile(res_sim[[3]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_2, col=col3_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,3], col=magma(20)[14], lwd=2)

hist3_3<-hist(res_sim[[3]][3,], breaks=seq(from=min(res_sim[[3]][3,])-0.1, 
                                           to=max(res_sim[[3]][3,])+0.1, by=0.008), plot=F)
col3_3<-ifelse(hist3_3$breaks < quantile(res_sim[[3]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_3$breaks >= quantile(res_sim[[3]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_3, col=col3_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,3], col=magma(20)[14], lwd=2)

hist3_4<-hist(res_sim[[3]][4,], breaks=seq(from=min(res_sim[[3]][4,])-0.1, 
                                           to=max(res_sim[[3]][4,])+0.2, by=0.01), plot=F)
col<-ifelse(hist3_4$breaks < quantile(res_sim[[3]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_4$breaks >= quantile(res_sim[[3]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,3], col=magma(20)[14], lwd=2)

hist3_5<-hist(res_sim[[3]][5,], breaks=seq(from=min(res_sim[[3]][5,])-0.3, 
                                           to=max(res_sim[[3]][5,])+0.1, by=0.01), plot=F)
col3_5<-ifelse(hist3_5$breaks < quantile(res_sim[[3]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist3_5$breaks >= quantile(res_sim[[3]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist3_5, col=col3_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,3], col=magma(20)[14], lwd=2)



hist4_1<-hist(res_sim[[4]][1,], breaks=seq(from=min(res_sim[[4]][1,])-0.00001, 
                                           to=max(res_sim[[4]][1,])+0.00001, by=0.0000015), plot=F)
col4_1<-ifelse(hist4_1$breaks < quantile(res_sim[[4]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist4_1$breaks >= 
                                               quantile(res_sim[[4]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist4_1, col=col4_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,4], col=magma(20)[14], lwd=2)

hist4_2<-hist(res_sim[[4]][2,], breaks=seq(from=min(res_sim[[4]][2,])-0.1,
                                           to=max(res_sim[[4]][2,])+0.1, by=0.01), plot=F)
col4_2<-ifelse(hist4_2$breaks < quantile(res_sim[[4]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_2$breaks >= quantile(res_sim[[4]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_2, col=col4_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,4], col=magma(20)[14], lwd=2)

hist4_3<-hist(res_sim[[4]][3,], breaks=seq(from=min(res_sim[[4]][3,])-0.1, 
                                           to=max(res_sim[[4]][3,])+0.1, by=0.008), plot=F)
col4_3<-ifelse(hist4_3$breaks < quantile(res_sim[[4]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_3$breaks >= quantile(res_sim[[4]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_3, col=col4_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,4], col=magma(20)[14], lwd=2)

hist4_4<-hist(res_sim[[4]][4,], breaks=seq(from=min(res_sim[[4]][4,])-0.1, 
                                           to=max(res_sim[[4]][4,])+0.2, by=0.01), plot=F)
col<-ifelse(hist4_4$breaks < quantile(res_sim[[4]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_4$breaks >= quantile(res_sim[[4]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,4], col=magma(20)[14], lwd=2)

hist4_5<-hist(res_sim[[4]][5,], breaks=seq(from=min(res_sim[[4]][5,])-0.4, 
                                           to=max(res_sim[[4]][5,])+0.1, by=0.01), plot=F)
col4_5<-ifelse(hist4_5$breaks < quantile(res_sim[[4]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist4_5$breaks >= quantile(res_sim[[4]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist4_5, col=col4_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,4], col=magma(20)[14], lwd=2)



hist5_1<-hist(res_sim[[5]][1,], breaks=seq(from=min(res_sim[[5]][1,])-0.00001, 
                                           to=max(res_sim[[5]][1,])+0.00001, by=0.0000025), plot=F)
col5_1<-ifelse(hist5_1$breaks < quantile(res_sim[[5]][1,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1), ifelse (hist5_1$breaks >= 
                                               quantile(res_sim[[5]][1,], probs=0.975), 
                                             rgb(0.1,0.1,0.1,0.1) , rgb(0.5,0.5,0.5,0.5) ))
plot(hist5_1, col=col5_1, border=F, main="", xlab=expression(paste("log(", rho, ")")))
abline(v=birds[1,5], col=magma(20)[14], lwd=2)

hist5_2<-hist(res_sim[[5]][2,], breaks=seq(from=min(res_sim[[5]][2,])-0.1,
                                           to=max(res_sim[[5]][2,])+0.1, by=0.01), plot=F)
col5_2<-ifelse(hist5_2$breaks < quantile(res_sim[[5]][2,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_2$breaks >= quantile(res_sim[[5]][2,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_2, col=col5_2, border=F, main="", xlab="mean(log(range size))")
abline(v=birds[2,5], col=magma(20)[14], lwd=2)

hist5_3<-hist(res_sim[[5]][3,], breaks=seq(from=min(res_sim[[5]][3,])-0.1, 
                                           to=max(res_sim[[5]][3,])+0.1, by=0.008), plot=F)
col5_3<-ifelse(hist5_3$breaks < quantile(res_sim[[5]][3,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_3$breaks >= quantile(res_sim[[5]][3,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_3, col=col5_3, border=F, main="", xlab="variance(log(range size))")
abline(v=birds[3,5], col=magma(20)[14], lwd=2)

hist5_4<-hist(res_sim[[5]][4,], breaks=seq(from=min(res_sim[[5]][4,])-0.1, 
                                           to=max(res_sim[[5]][4,])+0.2, by=0.02), plot=F)
col<-ifelse(hist5_4$breaks < quantile(res_sim[[5]][4,], probs=0.025), 
            rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_4$breaks >= quantile(res_sim[[5]][4,], probs=0.975), 
                                           rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_4, col=col, border=F, main="", xlab="mean(log(range overlap))")
abline(v=birds[4,5], col=magma(20)[14], lwd=2)

hist5_5<-hist(res_sim[[5]][5,], breaks=seq(from=min(res_sim[[5]][5,])-0.4, 
                                           to=max(res_sim[[5]][5,])+0.1, by=0.02), plot=F)
col5_5<-ifelse(hist5_5$breaks < quantile(res_sim[[5]][5,], probs=0.025), 
               rgb(0.1,0.1,0.1,0.1) , ifelse (hist5_5$breaks >= quantile(res_sim[[5]][5,], probs=0.975), 
                                              rgb(0.1,0.1,0.1,0.1) , rgb(0.6,0.6,0.6,0.6) ))
plot(hist5_5, col=col5_5, border=F, main="", xlab="variance(log(range overlap))")
abline(v=birds[5,5], col=magma(20)[14], lwd=2)

dev.off()



net.birds<-list()
net.birds$MAH<-matrix(nrow=S_MAH, ncol=3)
net.birds$CDH<-matrix(nrow=S_CDH, ncol=3)
net.birds$CAN<-matrix(nrow=S_CAN, ncol=3)
net.birds$NAN<-matrix(nrow=S_NAN, ncol=3)
net.birds$TEP<-matrix(nrow=S_TEP, ncol=3)

overlap.MAH[overlap.MAH > 0]<-overlap.CDH[overlap.CDH > 0]<-
  overlap.CAN[overlap.CAN > 0]<-overlap.NAN[overlap.NAN > 0]<-
  overlap.TEP[overlap.TEP > 0]<-1

net.birds$MAH[,1]<-sort(degree(overlap.MAH, gmode="graph"), decreasing=T)
net.birds$MAH[,2]<-sort(betweenness(overlap.MAH, gmode="graph"), decreasing=T)
net.birds$MAH[,3]<-sort(closeness(overlap.MAH, gmode="graph"), decreasing=T)

net.birds$CDH[,1]<-sort(degree(overlap.CDH, gmode="graph"), decreasing=T)
net.birds$CDH[,2]<-sort(betweenness(overlap.CDH, gmode="graph"), decreasing=T)
net.birds$CDH[,3]<-sort(closeness(overlap.CDH, gmode="graph"), decreasing=T)

net.birds$CAN[,1]<-sort(degree(overlap.CAN, gmode="graph"), decreasing=T)
net.birds$CAN[,2]<-sort(betweenness(overlap.CAN, gmode="graph"), decreasing=T)
net.birds$CAN[,3]<-sort(closeness(overlap.CAN, gmode="graph"), decreasing=T)

net.birds$NAN[,1]<-sort(degree(overlap.NAN, gmode="graph"), decreasing=T)
net.birds$NAN[,2]<-sort(betweenness(overlap.NAN, gmode="graph"), decreasing=T)
net.birds$NAN[,3]<-sort(closeness(overlap.NAN, gmode="graph"), decreasing=T)

net.birds$TEP[,1]<-sort(degree(overlap.TEP, gmode="graph"), decreasing=T)
net.birds$TEP[,2]<-sort(betweenness(overlap.TEP, gmode="graph"), decreasing=T)
net.birds$TEP[,3]<-sort(closeness(overlap.TEP, gmode="graph"), decreasing=T)

net.sim<-list()
for(i in 1:100){
  print(i)
  net.sim[[i]]<-list()
  
  net.sim[[i]]$MAH<-matrix(nrow=S_MAH, ncol=3)
  net.sim[[i]]$CDH<-matrix(nrow=S_CDH, ncol=3)
  net.sim[[i]]$CAN<-matrix(nrow=S_CAN, ncol=3)
  net.sim[[i]]$NAN<-matrix(nrow=S_NAN, ncol=3)
  net.sim[[i]]$TEP<-matrix(nrow=S_TEP, ncol=3)
  
  sim.MAH[[i]]$over[sim.MAH[[i]]$over > 0]<-1
  sim.CDH[[i]]$over[sim.CDH[[i]]$over > 0]<-1
  sim.CAN[[i]]$over[sim.CAN[[i]]$over > 0]<-1
  sim.NAN[[i]]$over[sim.NAN[[i]]$over > 0]<-1
  sim.TEP[[i]]$over[sim.TEP[[i]]$over > 0]<-1
  
  net.sim[[i]]$MAH[,1]<-sort(degree(sim.MAH[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$MAH[,2]<-sort(betweenness(sim.MAH[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$MAH[,3]<-sort(closeness(sim.MAH[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$CDH[,1]<-sort(degree(sim.CDH[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CDH[,2]<-sort(betweenness(sim.CDH[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CDH[,3]<-sort(closeness(sim.CDH[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$CAN[,1]<-sort(degree(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CAN[,2]<-sort(betweenness(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$CAN[,3]<-sort(closeness(sim.CAN[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$NAN[,1]<-sort(degree(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$NAN[,2]<-sort(betweenness(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$NAN[,3]<-sort(closeness(sim.NAN[[i]]$over, gmode="graph"), decreasing=T)

  net.sim[[i]]$TEP[,1]<-sort(degree(sim.TEP[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$TEP[,2]<-sort(betweenness(sim.TEP[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$TEP[,3]<-sort(closeness(sim.TEP[[i]]$over, gmode="graph"), decreasing=T)
}

load(file="data/net.sim.NT.RData")
pdf("figures/Figure7A.pdf")
layout(matrix(1:15, ncol=3, byrow=TRUE))
par(oma=c(2, 2, 2, 2), mar=c(4, 4, 0.5, 0.5))

rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col <- rbPal1(5)[as.numeric(cut(1:5, breaks = 5))]

plot(net.sim[[1]]$MAH[,1], main="", ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,1], col="gray")
}
lines(net.birds$MAH[,1], col=col[1]) 

plot(net.sim[[1]]$MAH[,2], main="", ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,2], col="gray")
}
lines(net.birds$MAH[,2], col=col[1]) 

plot(net.sim[[1]]$MAH[,3], main="", ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,3], col="gray")
}
lines(net.birds$MAH[,3], col=col[1]) 


plot(net.sim[[1]]$CDH[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,1], col="gray")
}
lines(net.birds$CDH[,1], col=col[2]) 

plot(net.sim[[1]]$CDH[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,2], col="gray")
}
lines(net.birds$CDH[,2], col=col[2]) 

plot(net.sim[[1]]$CDH[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,3], col="gray")
}
lines(net.birds$CDH[,3], col=col[2]) 


plot(net.sim[[1]]$CAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$CAN[,1], col="gray")
}
lines(net.birds$CAN[,1], col=col[3]) 

plot(net.sim[[1]]$CAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CAN[,2], col="gray")
}
lines(net.birds$CAN[,2], col=col[3]) 

plot(net.sim[[1]]$CAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CAN[,3], col="gray")
}
lines(net.birds$CAN[,3], col=col[3]) 


plot(net.sim[[1]]$NAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,1], col="gray")
}
lines(net.birds$NAN[,1], col=col[4]) 

plot(net.sim[[1]]$NAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,2], col="gray")
}
lines(net.birds$NAN[,2], col=col[4]) 

plot(net.sim[[1]]$NAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,3], col="gray")
}
lines(net.birds$NAN[,3], col=col[4]) 


plot(net.sim[[1]]$TEP[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,1], col="gray")
}
lines(net.birds$TEP[,1], col=col[5]) 

plot(net.sim[[1]]$TEP[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,2], col="gray")
}
lines(net.birds$TEP[,2], col=col[5]) 

plot(net.sim[[1]]$TEP[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,3], col="gray")
}
lines(net.birds$TEP[,3], col=col[5]) 

dev.off()

load(file="data/net.sim.T.RData")
pdf("figures/Figure7B.pdf")
layout(matrix(1:15, ncol=3, byrow=TRUE))
par(oma=c(2, 2, 2, 2), mar=c(4, 4, 0.5, 0.5))

rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col <- rbPal1(5)[as.numeric(cut(1:5, breaks = 5))]

plot(net.sim[[1]]$MAH[,1], main="", ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,1], col="gray")
}
lines(net.birds$MAH[,1], col=col[1]) 

plot(net.sim[[1]]$MAH[,2], main="", ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,2], col="gray")
}
lines(net.birds$MAH[,2], col=col[1]) 

plot(net.sim[[1]]$MAH[,3], main="", ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$MAH[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$MAH[,3], col="gray")
}
lines(net.birds$MAH[,3], col=col[1]) 


plot(net.sim[[1]]$CDH[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,1], col="gray")
}
lines(net.birds$CDH[,1], col=col[2]) 

plot(net.sim[[1]]$CDH[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,2], col="gray")
}
lines(net.birds$CDH[,2], col=col[2]) 

plot(net.sim[[1]]$CDH[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CDH[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CDH[,3], col="gray")
}
lines(net.birds$CDH[,3], col=col[2]) 


plot(net.sim[[1]]$CAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$CAN[,1], col="gray")
}
lines(net.birds$CAN[,1], col=col[3]) 

plot(net.sim[[1]]$CAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CAN[,2], col="gray")
}
lines(net.birds$CAN[,2], col=col[3]) 

plot(net.sim[[1]]$CAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$CAN[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$CAN[,3], col="gray")
}
lines(net.birds$CAN[,3], col=col[3]) 


plot(net.sim[[1]]$NAN[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,1], col="gray")
}
lines(net.birds$NAN[,1], col=col[4]) 

plot(net.sim[[1]]$NAN[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,2], col="gray")
}
lines(net.birds$NAN[,2], col=col[4]) 

plot(net.sim[[1]]$NAN[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$NAN[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$NAN[,3], col="gray")
}
lines(net.birds$NAN[,3], col=col[4]) 


plot(net.sim[[1]]$TEP[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,1]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,1], col="gray")
}
lines(net.birds$TEP[,1], col=col[5]) 

plot(net.sim[[1]]$TEP[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,2]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,2], col="gray")
}
lines(net.birds$TEP[,2], col=col[5]) 

plot(net.sim[[1]]$TEP[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$TEP[,3]), col="gray") 
for(i in 2:100){
  lines(net.sim[[i]]$TEP[,3], col="gray")
}
lines(net.birds$TEP[,3], col=col[5]) 

dev.off()