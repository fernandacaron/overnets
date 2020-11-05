rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

library(sna)
library(viridis)

source("codes/functions.R")

dat_a<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")[,c(10,12)]
rownames(dat_a)<-paste(
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")[,3],
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")[,4], sep="_")
dat_a[,c(1:2)]<-apply(dat_a[,c(1:2)], 2, function(x) as.numeric((x)))

dat_c<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/cdata.csv")[,c(10,12)]
rownames(dat_c)<-paste(
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/cdata.csv")[,3],
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/cdata.csv")[,4], sep="_")
dat_c[,c(1:2)]<-apply(dat_c[,c(1:2)], 2, function(x) as.numeric((x)))

dat_d<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/ddata.csv")[,c(8,10)]
rownames(dat_d)<-paste(
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/ddata.csv")[,3],
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/ddata.csv")[,4], sep="_")
dat_d[,c(1:2)]<-apply(dat_d[,c(1:2)], 2, function(x) as.numeric((x)))

dat_e<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/edata.csv")[,c(10,12)]
rownames(dat_e)<-paste(
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/edata.csv")[,1],
  read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/edata.csv")[,2], sep="_")
dat_e[,c(1:2)]<-apply(dat_e[,c(1:2)], 2, function(x) as.numeric((x)))

dat_a<-na.omit(dat_a)
dat_c<-na.omit(dat_c)
dat_d<-na.omit(dat_d)
dat_e<-na.omit(dat_e)

Sa<-nrow(dat_a)
Sc<-nrow(dat_c)
Sd<-nrow(dat_d)
Se<-nrow(dat_e)

range.sizes.a<-numeric()
for(i in 1:Sa){
  range.sizes.a<-dat_a$MAX-dat_a$MIN
}
names(range.sizes.a)<-rownames(dat_a)

range.sizes.c<-numeric()
for(i in 1:Sc){
  range.sizes.c<-dat_c$MAX-dat_c$MIN
}
names(range.sizes.c)<-rownames(dat_c)

range.sizes.d<-numeric()
for(i in 1:Sd){
  range.sizes.d<-dat_d$MAX-dat_d$MIN
}
names(range.sizes.d)<-rownames(dat_d)

range.sizes.e<-numeric()
for(i in 1:Se){
  range.sizes.e<-dat_e$MAX-dat_e$MIN
}
names(range.sizes.e)<-rownames(dat_e)

overlap.a<-matrix(nrow=Sa, ncol=Sa)
for (i in 1:Sa){
  for (j in 1:Sa){
    pair<-rbind(dat_a[i, ], dat_a[j, ])
    overlap.a[i, j]<-overcalc(pair[1, ], pair[2, ])
  }
}

overlap.c<-matrix(nrow=Sc, ncol=Sc)
for (i in 1:Sc){
  for (j in 1:Sc){
    pair<-rbind(dat_c[i, ], dat_c[j, ])
    overlap.c[i, j]<-overcalc(pair[1, ], pair[2, ])
  }
}

overlap.d<-matrix(nrow=Sd, ncol=Sd)
for (i in 1:Sd){
  for (j in 1:Sd){
    pair<-rbind(dat_d[i, ], dat_d[j, ])
    overlap.d[i, j]<-overcalc(pair[1, ], pair[2, ])
  }
}

overlap.e<-matrix(nrow=Se, ncol=Se)
for (i in 1:Se){
  for (j in 1:Se){
    pair<-rbind(dat_e[i, ], dat_e[j, ])
    overlap.e[i, j]<-overcalc(pair[1, ], pair[2, ])
  }
}

birds<-matrix(nrow=5, ncol=4)
colnames(birds)<-c("a", "c", "d", "e")
rownames(birds)<-c("rho", "log.mean.range.size",
                   "log.variance.range.size", 
                   "log.mean.range.overlap",
                   "log.variance.range.overlap")


birds[1,1]<-mean(range.sizes.a)/(max(dat_a$MAX)-min(dat_a$MIN))
birds[1,2]<-mean(range.sizes.c)/(max(dat_c$MAX)-min(dat_c$MIN))
birds[1,3]<-mean(range.sizes.d)/(max(dat_d$MAX)-min(dat_d$MIN))
birds[1,4]<-mean(range.sizes.e)/(max(dat_e$MAX)-min(dat_e$MIN))

birds[2,1]<-mean(log(diag(overlap.a)))
birds[2,2]<-mean(log(diag(overlap.c)))
birds[2,3]<-mean(log(diag(overlap.d)))
birds[2,4]<-mean(log(diag(overlap.e)))

birds[3,1]<-var(log(diag(overlap.a)))
birds[3,2]<-var(log(diag(overlap.c)))
birds[3,3]<-var(log(diag(overlap.d)))
birds[3,4]<-var(log(diag(overlap.e)))

diag(overlap.a)<-diag(overlap.c)<-diag(overlap.d)<-diag(overlap.e)<-NA

over.a<-overlap.a[overlap.a > 0]
over.c<-overlap.c[overlap.c > 0]
over.d<-overlap.d[overlap.d > 0]
over.e<-overlap.e[overlap.e > 0]

birds[4,1]<-mean(log(over.a), na.rm = TRUE)
birds[4,2]<-mean(log(over.c), na.rm = TRUE)
birds[4,3]<-mean(log(over.d), na.rm = TRUE)
birds[4,4]<-mean(log(over.e), na.rm = TRUE)

birds[5,1]<-var(as.vector(log(over.a)), na.rm = TRUE)
birds[5,2]<-var(as.vector(log(over.c)), na.rm = TRUE)
birds[5,3]<-var(as.vector(log(over.d)), na.rm = TRUE)
birds[5,4]<-var(as.vector(log(over.e)), na.rm = TRUE)

meanlog.a<-mean(log(range.sizes.a))
meanlog.c<-mean(log(range.sizes.c))
meanlog.d<-mean(log(range.sizes.d))
meanlog.e<-mean(log(range.sizes.e))

sim.a<-sim.c<-sim.d<-sim.e<-list()
for(i in 1:1000){
  sim.a[[i]]<-simRangeOver2(meanlog.a, Sa, truncate=TRUE, D_min=min(dat_a$MIN),
                           D_max=max(dat_a$MAX))
  sim.c[[i]]<-simRangeOver2(meanlog.c, Sc, truncate=TRUE, D_min=min(dat_c$MIN),
                          D_max=max(dat_c$MAX))
  sim.d[[i]]<-simRangeOver2(meanlog.d, Sd, truncate=TRUE, D_min=min(dat_d$MIN),
                          D_max=max(dat_d$MAX))
  sim.e[[i]]<-simRangeOver2(meanlog.e, Se, truncate=TRUE, D_min=min(dat_e$MIN),
                          D_max=max(dat_e$MAX))
}

res_sim<-list()
res_sim[[1]]<-res_sim[[2]]<-res_sim[[3]]<-res_sim[[4]]<-matrix(nrow=5, ncol=length(sim.a))
rownames(res_sim[[1]])<-rownames(res_sim[[2]])<-rownames(res_sim[[3]])<-
  rownames(res_sim[[4]])<-c("rho", "log.mean.range.size", 
                            "log.variance.range.size", 
                            "log.mean.range.overlap", 
                            "log.variance.range.overlap")
sim.a.over<-sim.c.over<-sim.d.over<-sim.e.over<-list()
for(i in 1:length(sim.a)){
  res_sim[[1]][1,i]<-sim.a[[i]]$rho
  res_sim[[2]][1,i]<-sim.c[[i]]$rho
  res_sim[[3]][1,i]<-sim.d[[i]]$rho
  res_sim[[4]][1,i]<-sim.e[[i]]$rho

  res_sim[[1]][2,i]<-mean(log(diag(sim.a[[i]]$over)))
  res_sim[[2]][2,i]<-mean(log(diag(sim.c[[i]]$over)))
  res_sim[[3]][2,i]<-mean(log(diag(sim.d[[i]]$over)))
  res_sim[[4]][2,i]<-mean(log(diag(sim.e[[i]]$over)))

  res_sim[[1]][3,i]<-var(log(diag(sim.a[[i]]$over)))
  res_sim[[2]][3,i]<-var(log(diag(sim.c[[i]]$over)))
  res_sim[[3]][3,i]<-var(log(diag(sim.d[[i]]$over)))
  res_sim[[4]][3,i]<-var(log(diag(sim.e[[i]]$over)))

  diag(sim.a[[i]]$over)<-diag(sim.c[[i]]$over)<-diag(sim.d[[i]]$over)<-
    diag(sim.e[[i]]$over)<-NA
  
  sim.a.over[[i]]<-sim.a[[i]]$over[sim.a[[i]]$over > 0]
  sim.c.over[[i]]<-sim.c[[i]]$over[sim.c[[i]]$over > 0]
  sim.d.over[[i]]<-sim.d[[i]]$over[sim.d[[i]]$over > 0]
  sim.e.over[[i]]<-sim.e[[i]]$over[sim.e[[i]]$over > 0]

  res_sim[[1]][4,i]<-mean(log(sim.a.over[[i]]), na.rm = TRUE)
  res_sim[[2]][4,i]<-mean(log(sim.c.over[[i]]), na.rm = TRUE)
  res_sim[[3]][4,i]<-mean(log(sim.d.over[[i]]), na.rm = TRUE)
  res_sim[[4]][4,i]<-mean(log(sim.e.over[[i]]), na.rm = TRUE)

  res_sim[[1]][5,i]<-var(as.vector(log(sim.a.over[[i]])), na.rm = TRUE)
  res_sim[[2]][5,i]<-var(as.vector(log(sim.c.over[[i]])), na.rm = TRUE)
  res_sim[[3]][5,i]<-var(as.vector(log(sim.d.over[[i]])), na.rm = TRUE)
  res_sim[[4]][5,i]<-var(as.vector(log(sim.e.over[[i]])), na.rm = TRUE)
}

pdf("figures/oi2.pdf")
layout(matrix(1:25, ncol=5, byrow = T))
par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

hist(res_sim[[1]][1,], col="white", border="gray", xlab="", main=expression(rho), breaks=seq(from=min(res_sim[[1]][1,])-0.1, to=max(res_sim[[1]][1,])+0.1, by=0.02))
sim1_rho_ci<-res_sim[[1]][1,][res_sim[[1]][1,] >= quantile(res_sim[[1]][1,], probs=0.025)]
sim1_rho_ci<-sim1_rho_ci[sim1_rho_ci <= quantile(res_sim[[1]][1,], probs=0.975)]
hist(sim1_rho_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][1,])-0.1, to=max(res_sim[[1]][1,])+0.1, by=0.02))
abline(v=birds[1,1], col="red")

hist(res_sim[[1]][2,], col="white", border="gray", xlab="", main="log(mean RS)", breaks=seq(from=min(res_sim[[1]][2,])-0.1, to=max(res_sim[[1]][2,])+0.1, by=0.02))
sim1_meanrs_ci<-res_sim[[1]][2,][res_sim[[1]][2,] >= quantile(res_sim[[1]][2,], probs=0.025)]
sim1_meanrs_ci<-sim1_meanrs_ci[sim1_meanrs_ci <= quantile(res_sim[[1]][2,], probs=0.975)]
hist(sim1_meanrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][2,])-0.1, to=max(res_sim[[1]][2,])+0.1, by=0.02))
abline(v=birds[2,1], col="red")

hist(res_sim[[1]][3,], col="white", border="gray", xlab="", main="log(variance RS)", breaks=seq(from=min(res_sim[[1]][3,])-0.1, to=max(res_sim[[1]][3,])+0.5, by=0.02))
sim1_varrs_ci<-res_sim[[1]][3,][res_sim[[1]][3,] >= quantile(res_sim[[1]][3,], probs=0.025)]
sim1_varrs_ci<-sim1_varrs_ci[sim1_varrs_ci<=quantile(res_sim[[1]][3,], probs=0.975)]
hist(sim1_varrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][3,])-0.1, to=max(res_sim[[1]][3,])+0.5, by=0.02))
abline(v=birds[3,1], col="red")

hist(res_sim[[1]][4,], col="white", border="gray", xlab="", main="log(mean RO)", breaks=seq(from=min(res_sim[[1]][4,])-0.1, to=max(res_sim[[1]][4,])+0.1, by=0.02))
sim1_meanro_ci<-res_sim[[1]][4,][res_sim[[1]][4,] >= quantile(res_sim[[1]][4,], probs=0.025)]
sim1_meanro_ci<-sim1_meanro_ci[sim1_meanro_ci <= quantile(res_sim[[1]][4,], probs=0.975)]
hist(sim1_meanro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][4,])-0.1, to=max(res_sim[[1]][4,])+0.1, by=0.02))
abline(v=birds[4,1], col="red")

hist(res_sim[[1]][5,], col="white", border="gray", xlab="", main="log(variance RO)", breaks=seq(from=min(res_sim[[1]][5,])-0.1, to=max(res_sim[[1]][5,])+1, by=0.02))
sim1_varro_ci<-res_sim[[1]][5,][res_sim[[1]][5,] >= quantile(res_sim[[1]][5,], probs=0.025)]
sim1_varro_ci<-sim1_varro_ci[sim1_varro_ci <= quantile(res_sim[[1]][5,], probs=0.975)]
hist(sim1_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[1]][5,])-0.1, to=max(res_sim[[1]][5,])+1, by=0.02))
abline(v=birds[5,1], col="red")



hist(res_sim[[2]][1,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[2]][1,])-0.1, to=max(res_sim[[2]][1,])+0.1, by=0.02))
sim2_rho_ci<-res_sim[[2]][1,][res_sim[[2]][1,] >= quantile(res_sim[[2]][1,], probs=0.025)]
sim2_rho_ci<-sim2_rho_ci[sim2_rho_ci <= quantile(res_sim[[2]][1,], probs=0.975)]
hist(sim2_rho_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[2]][1,])-0.1, to=max(res_sim[[2]][1,])+0.1, by=0.02))
abline(v=birds[1,2], col="red")

hist(res_sim[[2]][2,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[2]][2,])-0.1, to=max(res_sim[[2]][2,])+0.1, by=0.02))
sim2_meanrs_ci<-res_sim[[2]][2,][res_sim[[2]][2,] >= quantile(res_sim[[2]][2,], probs=0.025)]
sim2_meanrs_ci<-sim2_meanrs_ci[sim2_meanrs_ci <= quantile(res_sim[[2]][2,], probs=0.975)]
hist(sim2_meanrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[2]][2,])-0.1, to=max(res_sim[[2]][2,])+0.1, by=0.02))
abline(v=birds[2,2], col="red")

hist(res_sim[[2]][3,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[2]][3,])-0.1, to=max(res_sim[[2]][3,])+0.1, by=0.02))
sim2_varrs_ci<-res_sim[[2]][3,][res_sim[[2]][3,] >= quantile(res_sim[[2]][3,], probs=0.025)]
sim2_varrs_ci<-sim2_varrs_ci[sim2_varrs_ci<=quantile(res_sim[[2]][3,], probs=0.975)]
hist(sim2_varrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[2]][3,])-0.1, to=max(res_sim[[2]][3,])+0.1, by=0.02))
abline(v=birds[3,2], col="red")

hist(res_sim[[2]][4,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[2]][4,])-0.1, to=max(res_sim[[2]][4,])+0.1, by=0.02))
sim2_meanro_ci<-res_sim[[2]][4,][res_sim[[2]][4,] >= quantile(res_sim[[2]][4,], probs=0.025)]
sim2_meanro_ci<-sim2_meanro_ci[sim2_meanro_ci <= quantile(res_sim[[2]][4,], probs=0.975)]
hist(sim2_meanro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[2]][4,])-0.1, to=max(res_sim[[2]][4,])+0.1, by=0.02))
abline(v=birds[4,2], col="red")

hist(res_sim[[2]][5,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[2]][5,])-0.1, to=max(res_sim[[2]][5,])+0.2, by=0.02))
sim2_varro_ci<-res_sim[[2]][5,][res_sim[[2]][5,] >= quantile(res_sim[[2]][5,], probs=0.025)]
sim2_varro_ci<-sim2_varro_ci[sim2_varro_ci <= quantile(res_sim[[2]][5,], probs=0.975)]
hist(sim2_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[2]][5,])-0.1, to=max(res_sim[[2]][5,])+0.2, by=0.02))
abline(v=birds[5,2], col="red")



hist(res_sim[[3]][1,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[3]][1,])-0.1, to=max(res_sim[[3]][1,])+0.1, by=0.02))
sim3_rho_ci<-res_sim[[3]][1,][res_sim[[3]][1,] >= quantile(res_sim[[3]][1,], probs=0.025)]
sim3_rho_ci<-sim3_rho_ci[sim3_rho_ci <= quantile(res_sim[[3]][1,], probs=0.975)]
hist(sim3_rho_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[3]][1,])-0.1, to=max(res_sim[[3]][1,])+0.1, by=0.02))
abline(v=birds[1,3], col="red")

hist(res_sim[[3]][2,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[3]][2,])-0.1, to=max(res_sim[[3]][2,])+0.1, by=0.02))
sim3_meanrs_ci<-res_sim[[3]][2,][res_sim[[3]][2,] >= quantile(res_sim[[3]][2,], probs=0.025)]
sim3_meanrs_ci<-sim3_meanrs_ci[sim3_meanrs_ci <= quantile(res_sim[[3]][2,], probs=0.975)]
hist(sim3_meanrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[3]][2,])-0.1, to=max(res_sim[[3]][2,])+0.1, by=0.02))
abline(v=birds[2,3], col="red")

hist(res_sim[[3]][3,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[3]][3,])-0.1, to=max(res_sim[[3]][3,])+0.1, by=0.02))
sim3_varrs_ci<-res_sim[[3]][3,][res_sim[[3]][3,] >= quantile(res_sim[[3]][3,], probs=0.025)]
sim3_varrs_ci<-sim3_varrs_ci[sim3_varrs_ci<=quantile(res_sim[[3]][3,], probs=0.975)]
hist(sim3_varrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[3]][3,])-0.1, to=max(res_sim[[3]][3,])+0.1, by=0.02))
abline(v=birds[3,3], col="red")

hist(res_sim[[3]][4,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[3]][4,])-0.1, to=max(res_sim[[3]][4,])+0.1, by=0.02))
sim3_meanro_ci<-res_sim[[3]][4,][res_sim[[3]][4,] >= quantile(res_sim[[3]][4,], probs=0.025)]
sim3_meanro_ci<-sim3_meanro_ci[sim3_meanro_ci <= quantile(res_sim[[3]][4,], probs=0.975)]
hist(sim3_meanro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[3]][4,])-0.1, to=max(res_sim[[3]][4,])+0.1, by=0.02))
abline(v=birds[4,3], col="red")

hist(res_sim[[3]][5,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[3]][5,])-0.1, to=max(res_sim[[3]][5,])+0.1, by=0.02))
sim3_varro_ci<-res_sim[[3]][5,][res_sim[[3]][5,] >= quantile(res_sim[[3]][5,], probs=0.025)]
sim3_varro_ci<-sim3_varro_ci[sim3_varro_ci <= quantile(res_sim[[3]][5,], probs=0.975)]
hist(sim3_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[3]][5,])-0.1, to=max(res_sim[[3]][5,])+0.1, by=0.02))
abline(v=birds[5,3], col="red")



hist(res_sim[[4]][1,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[4]][1,])-0.1, to=max(res_sim[[4]][1,])+0.1, by=0.02))
sim4_rho_ci<-res_sim[[4]][1,][res_sim[[4]][1,] >= quantile(res_sim[[4]][1,], probs=0.025)]
sim4_rho_ci<-sim4_rho_ci[sim4_rho_ci <= quantile(res_sim[[4]][1,], probs=0.975)]
hist(sim4_rho_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[4]][1,])-0.1, to=max(res_sim[[4]][1,])+0.1, by=0.02))
abline(v=birds[1,4], col="red")

hist(res_sim[[4]][2,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[4]][2,])-0.1, to=max(res_sim[[4]][2,])+0.1, by=0.02))
sim4_meanrs_ci<-res_sim[[4]][2,][res_sim[[4]][2,] >= quantile(res_sim[[4]][2,], probs=0.025)]
sim4_meanrs_ci<-sim4_meanrs_ci[sim4_meanrs_ci <= quantile(res_sim[[4]][2,], probs=0.975)]
hist(sim4_meanrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[4]][2,])-0.1, to=max(res_sim[[4]][2,])+0.1, by=0.02))
abline(v=birds[2,4], col="red")

hist(res_sim[[4]][3,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[4]][3,])-0.1, to=max(res_sim[[4]][3,])+0.1, by=0.02))
sim4_varrs_ci<-res_sim[[4]][3,][res_sim[[4]][3,] >= quantile(res_sim[[4]][3,], probs=0.025)]
sim4_varrs_ci<-sim4_varrs_ci[sim4_varrs_ci<=quantile(res_sim[[4]][3,], probs=0.975)]
hist(sim4_varrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[4]][3,])-0.1, to=max(res_sim[[4]][3,])+0.1, by=0.02))
abline(v=birds[3,4], col="red")

hist(res_sim[[4]][4,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[4]][4,])-0.1, to=max(res_sim[[4]][4,])+0.1, by=0.02))
sim4_meanro_ci<-res_sim[[4]][4,][res_sim[[4]][4,] >= quantile(res_sim[[4]][4,], probs=0.025)]
sim4_meanro_ci<-sim4_meanro_ci[sim4_meanro_ci <= quantile(res_sim[[4]][4,], probs=0.975)]
hist(sim4_meanro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[4]][4,])-0.1, to=max(res_sim[[4]][4,])+0.1, by=0.02))
abline(v=birds[4,4], col="red")

hist(res_sim[[4]][5,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[4]][5,])-0.1, to=max(res_sim[[4]][5,])+0.5, by=0.02))
sim4_varro_ci<-res_sim[[4]][5,][res_sim[[4]][5,] >= quantile(res_sim[[4]][5,], probs=0.025)]
sim4_varro_ci<-sim4_varro_ci[sim4_varro_ci <= quantile(res_sim[[4]][5,], probs=0.975)]
hist(sim4_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[4]][5,])-0.1, to=max(res_sim[[4]][5,])+0.5, by=0.02))
abline(v=birds[5,4], col="red")


dev.off()


net.birds<-list()
net.birds$a<-matrix(nrow=Sa, ncol=3)
net.birds$c<-matrix(nrow=Sc, ncol=3)
net.birds$d<-matrix(nrow=Sd, ncol=3)
net.birds$e<-matrix(nrow=Se, ncol=3)

overlap.a[overlap.a > 0] <- overlap.c[overlap.c > 0] <- overlap.d[overlap.d > 0] <- 
  overlap.e[overlap.e > 0] <- 1

net.birds$a[,1]<-sort(degree(overlap.a, gmode="graph"), decreasing=T)
net.birds$a[,2]<-sort(betweenness(overlap.a, gmode="graph"), decreasing=T)
net.birds$a[,3]<-sort(closeness(overlap.a, gmode="graph"), decreasing=T)

net.birds$c[,1]<-sort(degree(overlap.c, gmode="graph"), decreasing=T)
net.birds$c[,2]<-sort(betweenness(overlap.c, gmode="graph"), decreasing=T)
net.birds$c[,3]<-sort(closeness(overlap.c, gmode="graph"), decreasing=T)

net.birds$d[,1]<-sort(degree(overlap.d, gmode="graph"), decreasing=T)
net.birds$d[,2]<-sort(betweenness(overlap.d, gmode="graph"), decreasing=T)
net.birds$d[,3]<-sort(closeness(overlap.d, gmode="graph"), decreasing=T)

net.birds$e[,1]<-sort(degree(overlap.e, gmode="graph"), decreasing=T)
net.birds$e[,2]<-sort(betweenness(overlap.e, gmode="graph"), decreasing=T)
net.birds$e[,3]<-sort(closeness(overlap.e, gmode="graph"), decreasing=T)

net.sim<-list()
for(i in 1:1000){
  net.sim[[i]]<-list()
  
  net.sim[[i]]$a<-matrix(nrow=Sa, ncol=3)
  net.sim[[i]]$c<-matrix(nrow=Sc, ncol=3)
  net.sim[[i]]$d<-matrix(nrow=Sd, ncol=3)
  net.sim[[i]]$e<-matrix(nrow=Se, ncol=3)

  sim.a[[i]]$over[sim.a[[i]]$over > 0] <- sim.c[[i]]$over[sim.c[[i]]$over > 0] <- 
    sim.d[[i]]$over[sim.d[[i]]$over > 0] <- sim.e[[i]]$over[sim.e[[i]]$over > 0] <- 1
  
  net.sim[[i]]$a[,1]<-sort(degree(sim.a[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$a[,2]<-sort(betweenness(sim.a[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$a[,3]<-sort(closeness(sim.a[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$c[,1]<-sort(degree(sim.c[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$c[,2]<-sort(betweenness(sim.c[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$c[,3]<-sort(closeness(sim.c[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$d[,1]<-sort(degree(sim.d[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$d[,2]<-sort(betweenness(sim.d[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$d[,3]<-sort(closeness(sim.d[[i]]$over, gmode="graph"), decreasing=T)
  
  net.sim[[i]]$e[,1]<-sort(degree(sim.e[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$e[,2]<-sort(betweenness(sim.e[[i]]$over, gmode="graph"), decreasing=T)
  net.sim[[i]]$e[,3]<-sort(closeness(sim.e[[i]]$over, gmode="graph"), decreasing=T)
}

pdf("figures/oi.pdf")
layout(matrix(1:15, ncol=3, byrow=TRUE))
par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col <- rbPal1(4)[as.numeric(cut(1:4, breaks = 4))]

plot(net.sim[[1]]$a[,1], main="Degree", ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$a[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$a[,1], col="gray")
}
lines(net.birds$a[,1], col=col[1]) 

plot(net.sim[[1]]$a[,2], main="Betweenness", ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$a[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$a[,2], col="gray")
}
lines(net.birds$a[,2], col=col[1]) 

plot(net.sim[[1]]$a[,3], main="Closeness", ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$a[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$a[,3], col="gray")
}
lines(net.birds$a[,3], col=col[1]) 


plot(net.sim[[1]]$c[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$c[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$c[,1], col="gray")
}
lines(net.birds$c[,1], col=col[2]) 

plot(net.sim[[1]]$c[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$c[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$c[,2], col="gray")
}
lines(net.birds$c[,2], col=col[2]) 

plot(net.sim[[1]]$c[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$c[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$c[,3], col="gray")
}
lines(net.birds$c[,3], col=col[2]) 


plot(net.sim[[1]]$d[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$d[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$d[,1], col="gray")
}
lines(net.birds$d[,1], col=col[3]) 

plot(net.sim[[1]]$d[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$d[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$d[,2], col="gray")
}
lines(net.birds$d[,2], col=col[3]) 

plot(net.sim[[1]]$d[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$d[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$d[,3], col="gray")
}
lines(net.birds$d[,3], col=col[3]) 


plot(net.sim[[1]]$e[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$e[,1]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$e[,1], col="gray")
}
lines(net.birds$e[,1], col=col[4]) 

plot(net.sim[[1]]$e[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$e[,2]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$e[,2], col="gray")
}
lines(net.birds$e[,2], col=col[4]) 

plot(net.sim[[1]]$e[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$e[,3]), col="gray") 
for(i in 2:1000){
  lines(net.sim[[i]]$e[,3], col="gray")
}
lines(net.birds$e[,3], col=col[4]) 

dev.off()