rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

library(sna)
library(viridis)

source("codes/functions.R")

dat<-read.csv("data/birds_Hawkins.csv", row.names=1)
dat[,c(1:10)]<-apply(dat[,c(1:10)], 2, function(x) as.numeric((x)))

dat1<-dat[1:2]
dat2<-dat[3:4]
dat3<-dat[5:6]
dat4<-dat[7:8]
dat5<-dat[9:10]

dat1<-na.omit(dat1)
dat2<-na.omit(dat2)
dat3<-na.omit(dat3)
dat4<-na.omit(dat4)
dat5<-na.omit(dat5)

S1<-nrow(dat1)
S2<-nrow(dat2)
S3<-nrow(dat3)
S4<-nrow(dat4)
S5<-nrow(dat5)

range.sizes1<-numeric()
for(i in 1:S1){
  range.sizes1<-dat1$Anjanaharidesud_high-dat1$Anjanaharidesud_low
}
names(range.sizes1)<-rownames(dat1)

range.sizes2<-numeric()
for(i in 1:S2){
  range.sizes2<-dat2$Marojejy_high-dat2$Marojejy_low
}
names(range.sizes2)<-rownames(dat2)

range.sizes3<-numeric()
for(i in 1:S3){
  range.sizes3<-dat3$Zahamena_high-dat3$Zahamena_low
}
names(range.sizes3)<-rownames(dat3)

range.sizes4<-numeric()
for(i in 1:S4){
  range.sizes4<-dat4$Andringitra_high-dat4$Andringitra_low
}
names(range.sizes4)<-rownames(dat4)

range.sizes5<-numeric()
for(i in 1:S5){
  range.sizes5<-dat5$Andohahela_high-dat5$Andohahela_low
}
names(range.sizes5)<-rownames(dat5)

overlap1<-matrix(nrow=S1, ncol=S1)
for (i in 1:S1) {
  for (j in 1:S1) {
    pair <- rbind(dat1[i, ], dat1[j, ])
    overlap1[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

overlap2<-matrix(nrow=S2, ncol=S2)
for (i in 1:S2) {
  for (j in 1:S2) {
    pair <- rbind(dat2[i, ], dat2[j, ])
    overlap2[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

overlap3<-matrix(nrow=S3, ncol=S3)
for (i in 1:S3) {
  for (j in 1:S3) {
    pair <- rbind(dat3[i, ], dat3[j, ])
    overlap3[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

overlap4<-matrix(nrow=S4, ncol=S4)
for (i in 1:S4) {
  for (j in 1:S4) {
    pair <- rbind(dat4[i, ], dat4[j, ])
    overlap4[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

overlap5<-matrix(nrow=S5, ncol=S5)
for (i in 1:S5) {
  for (j in 1:S5) {
    pair <- rbind(dat5[i, ], dat5[j, ])
    overlap5[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

birds<-matrix(nrow=5, ncol=5)
colnames(birds)<-c("Anja", "Maro", "Zaha", "Andr", "Ando")
rownames(birds)<-c("rho", "log.mean.range.size",
                   "log.variance.range.size", 
                   "log.mean.range.overlap",
                   "log.variance.range.overlap")


birds[1,1]<-mean(range.sizes1)/(max(dat1$Anjanaharidesud_high)-min(dat1$Anjanaharidesud_low))
birds[1,2]<-mean(range.sizes2)/(max(dat2$Marojejy_high)-min(dat2$Marojejy_low))
birds[1,3]<-mean(range.sizes3)/(max(dat3$Zahamena_high)-min(dat3$Zahamena_low))
birds[1,4]<-mean(range.sizes4)/(max(dat4$Andringitra_high)-min(dat4$Andringitra_low))
birds[1,5]<-mean(range.sizes5)/(max(dat5$Andohahela_high)-min(dat5$Andohahela_low))

birds[2,1]<-mean(log(diag(overlap1)))
birds[2,2]<-mean(log(diag(overlap2)))
birds[2,3]<-mean(log(diag(overlap3)))
birds[2,4]<-mean(log(diag(overlap4)))
birds[2,5]<-mean(log(diag(overlap5)))

birds[3,1]<-var(log(diag(overlap1)))
birds[3,2]<-var(log(diag(overlap2)))
birds[3,3]<-var(log(diag(overlap3)))
birds[3,4]<-var(log(diag(overlap4)))
birds[3,5]<-var(log(diag(overlap5)))

diag(overlap1)<-diag(overlap2)<-diag(overlap3)<-diag(overlap4)<-diag(overlap5)<-NA

over1<-overlap1[overlap1 > 0]
over2<-overlap2[overlap2 > 0]
over3<-overlap3[overlap3 > 0]
over4<-overlap4[overlap4 > 0]
over5<-overlap5[overlap5 > 0]

birds[4,1]<-mean(log(over1), na.rm = TRUE)
birds[4,2]<-mean(log(over2), na.rm = TRUE)
birds[4,3]<-mean(log(over3), na.rm = TRUE)
birds[4,4]<-mean(log(over4), na.rm = TRUE)
birds[4,5]<-mean(log(over5), na.rm = TRUE)

birds[5,1]<-var(as.vector(log(over1)), na.rm = TRUE)
birds[5,2]<-var(as.vector(log(over2)), na.rm = TRUE)
birds[5,3]<-var(as.vector(log(over3)), na.rm = TRUE)
birds[5,4]<-var(as.vector(log(over4)), na.rm = TRUE)
birds[5,5]<-var(as.vector(log(over5)), na.rm = TRUE)

meanlog1<-mean(log(range.sizes1))
meanlog2<-mean(log(range.sizes2))
meanlog3<-mean(log(range.sizes3))
meanlog4<-mean(log(range.sizes4))
meanlog5<-mean(log(range.sizes5))

sim1<-sim2<-sim3<-sim4<-sim5<-list()
for(i in 1:1000){
  sim1[[i]]<-simRangeOver2(meanlog1, S1, truncate=TRUE, D_min=min(dat1$Anjanaharidesud_low),
                      D_max= max(dat1$Anjanaharidesud_high))
  sim2[[i]]<-simRangeOver2(meanlog2, S2, truncate=TRUE, D_min=min(dat2$Marojejy_low),
                      D_max= max(dat2$Marojejy_high))
  sim3[[i]]<-simRangeOver2(meanlog3, S3, truncate=TRUE, D_min=min(dat3$Zahamena_low),
                      D_max= max(dat3$Zahamena_high))
  sim4[[i]]<-simRangeOver2(meanlog4, S4, truncate=TRUE, D_min=min(dat4$Andringitra_low),
                      D_max= max(dat4$Andringitra_high))
  sim5[[i]]<-simRangeOver2(meanlog5, S5, truncate=TRUE, D_min=min(dat5$Andohahela_low),
                      D_max= max(dat5$Andohahela_high))
}

res_sim<-list()
res_sim[[1]]<-res_sim[[2]]<-res_sim[[3]]<-res_sim[[4]]<-res_sim[[5]]<-
  matrix(nrow=5, ncol=length(sim1))
rownames(res_sim[[1]])<-rownames(res_sim[[2]])<-rownames(res_sim[[3]])<-
  rownames(res_sim[[4]])<-rownames(res_sim[[5]])<-c("rho", "log.mean.range.size", 
                                                    "log.variance.range.size", 
                                                    "log.mean.range.overlap", 
                                                    "log.variance.range.overlap")
sim1.over<-sim2.over<-sim3.over<-sim4.over<-sim5.over<-list()
for(i in 1:length(sim1)){
  res_sim[[1]][1,i]<-sim1[[i]]$rho
  res_sim[[2]][1,i]<-sim2[[i]]$rho
  res_sim[[3]][1,i]<-sim3[[i]]$rho
  res_sim[[4]][1,i]<-sim4[[i]]$rho
  res_sim[[5]][1,i]<-sim5[[i]]$rho

  res_sim[[1]][2,i]<-mean(log(diag(sim1[[i]]$over)))
  res_sim[[2]][2,i]<-mean(log(diag(sim2[[i]]$over)))
  res_sim[[3]][2,i]<-mean(log(diag(sim3[[i]]$over)))
  res_sim[[4]][2,i]<-mean(log(diag(sim4[[i]]$over)))
  res_sim[[5]][2,i]<-mean(log(diag(sim5[[i]]$over)))

  res_sim[[1]][3,i]<-var(log(diag(sim1[[i]]$over)))
  res_sim[[2]][3,i]<-var(log(diag(sim2[[i]]$over)))
  res_sim[[3]][3,i]<-var(log(diag(sim3[[i]]$over)))
  res_sim[[4]][3,i]<-var(log(diag(sim4[[i]]$over)))
  res_sim[[5]][3,i]<-var(log(diag(sim5[[i]]$over)))

  diag(sim1[[i]]$over)<-diag(sim2[[i]]$over)<-diag(sim3[[i]]$over)<-
    diag(sim4[[i]]$over)<-diag(sim5[[i]]$over)<-NA

  sim1.over[[i]]<-sim1[[i]]$over[sim1[[i]]$over > 0]
  sim2.over[[i]]<-sim2[[i]]$over[sim2[[i]]$over > 0]
  sim3.over[[i]]<-sim3[[i]]$over[sim3[[i]]$over > 0]
  sim4.over[[i]]<-sim4[[i]]$over[sim4[[i]]$over > 0]
  sim5.over[[i]]<-sim5[[i]]$over[sim5[[i]]$over > 0]

  res_sim[[1]][4,i]<-mean(log(sim1.over[[i]]), na.rm = TRUE)
  res_sim[[2]][4,i]<-mean(log(sim2.over[[i]]), na.rm = TRUE)
  res_sim[[3]][4,i]<-mean(log(sim3.over[[i]]), na.rm = TRUE)
  res_sim[[4]][4,i]<-mean(log(sim4.over[[i]]), na.rm = TRUE)
  res_sim[[5]][4,i]<-mean(log(sim5.over[[i]]), na.rm = TRUE)

  res_sim[[1]][5,i]<-var(as.vector(log(sim1.over[[i]])), na.rm = TRUE)
  res_sim[[2]][5,i]<-var(as.vector(log(sim2.over[[i]])), na.rm = TRUE)
  res_sim[[3]][5,i]<-var(as.vector(log(sim3.over[[i]])), na.rm = TRUE)
  res_sim[[4]][5,i]<-var(as.vector(log(sim4.over[[i]])), na.rm = TRUE)
  res_sim[[5]][5,i]<-var(as.vector(log(sim5.over[[i]])), na.rm = TRUE)

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



hist(res_sim[[5]][1,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[5]][1,])-0.1, to=max(res_sim[[5]][1,])+0.1, by=0.02))
sim5_rho_ci<-res_sim[[5]][1,][res_sim[[5]][1,] >= quantile(res_sim[[5]][1,], probs=0.025)]
sim5_rho_ci<-sim5_rho_ci[sim5_rho_ci <= quantile(res_sim[[5]][1,], probs=0.975)]
hist(sim5_rho_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[5]][1,])-0.1, to=max(res_sim[[5]][1,])+0.1, by=0.02))
abline(v=birds[1,5], col="red")

hist(res_sim[[5]][2,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[5]][2,])-0.1, to=max(res_sim[[5]][2,])+0.1, by=0.02))
sim5_meanrs_ci<-res_sim[[5]][2,][res_sim[[5]][2,] >= quantile(res_sim[[5]][2,], probs=0.025)]
sim5_meanrs_ci<-sim5_meanrs_ci[sim5_meanrs_ci <= quantile(res_sim[[5]][2,], probs=0.975)]
hist(sim5_meanrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[5]][2,])-0.1, to=max(res_sim[[5]][2,])+0.1, by=0.02))
abline(v=birds[2,5], col="red")

hist(res_sim[[5]][3,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[5]][3,])-0.1, to=max(res_sim[[5]][3,])+0.1, by=0.02))
sim5_varrs_ci<-res_sim[[5]][3,][res_sim[[5]][3,] >= quantile(res_sim[[5]][3,], probs=0.025)]
sim5_varrs_ci<-sim5_varrs_ci[sim5_varrs_ci<=quantile(res_sim[[5]][3,], probs=0.975)]
hist(sim5_varrs_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[5]][3,])-0.1, to=max(res_sim[[5]][3,])+0.1, by=0.02))
abline(v=birds[3,5], col="red")

hist(res_sim[[5]][4,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[5]][4,])-0.1, to=max(res_sim[[5]][4,])+0.1, by=0.02))
sim5_meanro_ci<-res_sim[[5]][4,][res_sim[[5]][4,] >= quantile(res_sim[[5]][4,], probs=0.025)]
sim5_meanro_ci<-sim5_meanro_ci[sim5_meanro_ci <= quantile(res_sim[[5]][4,], probs=0.975)]
hist(sim5_meanro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[5]][4,])-0.1, to=max(res_sim[[5]][4,])+0.1, by=0.02))
abline(v=birds[4,5], col="red")

hist(res_sim[[5]][5,], col="white", border="gray", xlab="", main="", breaks=seq(from=min(res_sim[[5]][5,])-0.1, to=max(res_sim[[5]][5,])+0.2, by=0.02))
sim5_varro_ci<-res_sim[[5]][5,][res_sim[[5]][5,] >= quantile(res_sim[[5]][5,], probs=0.025)]
sim5_varro_ci<-sim5_varro_ci[sim5_varro_ci <= quantile(res_sim[[5]][5,], probs=0.975)]
hist(sim5_varro_ci, col="gray", border="gray", add=T, breaks=seq(from=min(res_sim[[5]][5,])-0.1, to=max(res_sim[[5]][5,])+0.2, by=0.02))
abline(v=birds[5,5], col="red")

dev.off()


net.birds<-list()
net.birds$Anja<-matrix(nrow=S1, ncol=3)
net.birds$Maro<-matrix(nrow=S2, ncol=3)
net.birds$Zaha<-matrix(nrow=S3, ncol=3)
net.birds$Andr<-matrix(nrow=S4, ncol=3)
net.birds$Ando<-matrix(nrow=S5, ncol=3)

overlap1[overlap1 > 0] <- overlap2[overlap2 > 0] <- overlap3[overlap3 > 0] <- 
  overlap4[overlap4 > 0] <- overlap5[overlap5 > 0] <- 1

net.birds$Anja[,1]<-sort(degree(overlap1, gmode="graph"), decreasing=T)
net.birds$Anja[,2]<-sort(betweenness(overlap1, gmode="graph"), decreasing=T)
net.birds$Anja[,3]<-sort(closeness(overlap1, gmode="graph"), decreasing=T)

net.birds$Maro[,1]<-sort(degree(overlap2, gmode="graph"), decreasing=T)
net.birds$Maro[,2]<-sort(betweenness(overlap2, gmode="graph"), decreasing=T)
net.birds$Maro[,3]<-sort(closeness(overlap2, gmode="graph"), decreasing=T)

net.birds$Zaha[,1]<-sort(degree(overlap3, gmode="graph"), decreasing=T)
net.birds$Zaha[,2]<-sort(betweenness(overlap3, gmode="graph"), decreasing=T)
net.birds$Zaha[,3]<-sort(closeness(overlap3, gmode="graph"), decreasing=T)

net.birds$Andr[,1]<-sort(degree(overlap4, gmode="graph"), decreasing=T)
net.birds$Andr[,2]<-sort(betweenness(overlap4, gmode="graph"), decreasing=T)
net.birds$Andr[,3]<-sort(closeness(overlap4, gmode="graph"), decreasing=T)

net.birds$Ando[,1]<-sort(degree(overlap5, gmode="graph"), decreasing=T)
net.birds$Ando[,2]<-sort(betweenness(overlap5, gmode="graph"), decreasing=T)
net.birds$Ando[,3]<-sort(closeness(overlap5, gmode="graph"), decreasing=T)

net.sim<-list()
for(i in 1:1000){
	net.sim[[i]]<-list()

	net.sim[[i]]$Anja<-matrix(nrow=S1, ncol=3)
	net.sim[[i]]$Maro<-matrix(nrow=S2, ncol=3)
	net.sim[[i]]$Zaha<-matrix(nrow=S3, ncol=3)
	net.sim[[i]]$Andr<-matrix(nrow=S4, ncol=3)
	net.sim[[i]]$Ando<-matrix(nrow=S5, ncol=3)

	sim1[[i]]$over[sim1[[i]]$over > 0] <- sim2[[i]]$over[sim2[[i]]$over > 0] <- sim3[[i]]$over[sim3[[i]]$over > 0] <- 
		sim4[[i]]$over[sim4[[i]]$over > 0] <- sim5[[i]]$over[sim5[[i]]$over > 0] <- 1

	net.sim[[i]]$Anja[,1]<-sort(degree(sim1[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Anja[,2]<-sort(betweenness(sim1[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Anja[,3]<-sort(closeness(sim1[[i]]$over, gmode="graph"), decreasing=T)

	net.sim[[i]]$Maro[,1]<-sort(degree(sim2[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Maro[,2]<-sort(betweenness(sim2[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Maro[,3]<-sort(closeness(sim2[[i]]$over, gmode="graph"), decreasing=T)

	net.sim[[i]]$Zaha[,1]<-sort(degree(sim3[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Zaha[,2]<-sort(betweenness(sim3[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Zaha[,3]<-sort(closeness(sim3[[i]]$over, gmode="graph"), decreasing=T)

	net.sim[[i]]$Andr[,1]<-sort(degree(sim4[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Andr[,2]<-sort(betweenness(sim4[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Andr[,3]<-sort(closeness(sim4[[i]]$over, gmode="graph"), decreasing=T)

	net.sim[[i]]$Ando[,1]<-sort(degree(sim5[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Ando[,2]<-sort(betweenness(sim5[[i]]$over, gmode="graph"), decreasing=T)
	net.sim[[i]]$Ando[,3]<-sort(closeness(sim5[[i]]$over, gmode="graph"), decreasing=T)
}

pdf("figures/oi.pdf")
layout(matrix(1:15, ncol=3, byrow=TRUE))
par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

rbPal1 <- colorRampPalette(c(viridis(20)[6:14]), alpha=TRUE)
col <- rbPal1(5)[as.numeric(cut(1:5, breaks = 5))]

plot(net.sim[[1]]$Anja[,1], main="Degree", ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$Anja[,1]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Anja[,1], col="gray")
}
lines(net.birds$Anja[,1], col=col[1]) 

plot(net.sim[[1]]$Anja[,2], main="Betweenness", ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Anja[,2]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Anja[,2], col="gray")
}
lines(net.birds$Anja[,2], col=col[1]) 

plot(net.sim[[1]]$Anja[,3], main="Closeness", ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Anja[,3]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Anja[,3], col="gray")
}
lines(net.birds$Anja[,3], col=col[1]) 


plot(net.sim[[1]]$Maro[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$Maro[,1]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Maro[,1], col="gray")
}
lines(net.birds$Maro[,1], col=col[2]) 

plot(net.sim[[1]]$Maro[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Maro[,2]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Maro[,2], col="gray")
}
lines(net.birds$Maro[,2], col=col[2]) 

plot(net.sim[[1]]$Maro[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Maro[,3]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Maro[,3], col="gray")
}
lines(net.birds$Maro[,3], col=col[2]) 


plot(net.sim[[1]]$Zaha[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$Zaha[,1]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Zaha[,1], col="gray")
}
lines(net.birds$Zaha[,1], col=col[3]) 

plot(net.sim[[1]]$Zaha[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Zaha[,2]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Zaha[,2], col="gray")
}
lines(net.birds$Zaha[,2], col=col[3]) 

plot(net.sim[[1]]$Zaha[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Zaha[,3]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Zaha[,3], col="gray")
}
lines(net.birds$Zaha[,3], col=col[3]) 


plot(net.sim[[1]]$Andr[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$Andr[,1]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Andr[,1], col="gray")
}
lines(net.birds$Andr[,1], col=col[4]) 

plot(net.sim[[1]]$Andr[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Andr[,2]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Andr[,2], col="gray")
}
lines(net.birds$Andr[,2], col=col[4]) 

plot(net.sim[[1]]$Andr[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Andr[,3]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Andr[,3], col="gray")
}
lines(net.birds$Andr[,3], col=col[4]) 

plot(net.sim[[1]]$Ando[,1], ylab="Degree", xlab="rank", type="l", ylim=range(net.sim[[1]]$Ando[,1]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Ando[,1], col="gray")
}
lines(net.birds$Ando[,1], col=col[5]) 

plot(net.sim[[1]]$Ando[,2], ylab="Betweenness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Ando[,2]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Ando[,2], col="gray")
}
lines(net.birds$Ando[,2], col=col[5]) 

plot(net.sim[[1]]$Ando[,3], ylab="Closeness", xlab="rank", type="l", ylim=range(net.sim[[1]]$Ando[,3]), col="gray") 
for(i in 2:1000){
	lines(net.sim[[i]]$Ando[,3], col="gray")
}
lines(net.birds$Ando[,3], col=col[5]) 

dev.off()