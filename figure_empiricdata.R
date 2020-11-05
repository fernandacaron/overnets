rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

library(sna)
library(igraph)

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

rho1<-mean(range.sizes1)/(max(dat1$Anjanaharidesud_high)-min(dat1$Anjanaharidesud_low))
rho2<-mean(range.sizes2)/(max(dat2$Marojejy_high)-min(dat2$Marojejy_low))
rho3<-mean(range.sizes3)/(max(dat3$Zahamena_high)-min(dat3$Zahamena_low))
rho4<-mean(range.sizes4)/(max(dat4$Andringitra_high)-min(dat4$Andringitra_low))
rho5<-mean(range.sizes5)/(max(dat5$Andohahela_high)-min(dat5$Andohahela_low))

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

network1<-graph_from_adjacency_matrix(overlap1, mode="undirected", weighted = T,
                                      diag=F)
network2<-graph_from_adjacency_matrix(overlap2, mode="undirected", weighted = T,
                                      diag=F)
network3<-graph_from_adjacency_matrix(overlap3, mode="undirected", weighted = T,
                                      diag=F)
network4<-graph_from_adjacency_matrix(overlap4, mode="undirected", weighted = T,
                                      diag=F)
network5<-graph_from_adjacency_matrix(overlap5, mode="undirected", weighted = T,
                                      diag=F)

l1<-layout.circle(network1)
l1<-layout.norm(l1, ymin=-1, ymax=1, xmin=-1, xmax=1)
l2<-layout.circle(network2)
l2<-layout.norm(l2, ymin=-1, ymax=1, xmin=-1, xmax=1)
l3<-layout.circle(network3)
l3<-layout.norm(l3, ymin=-1, ymax=1, xmin=-1, xmax=1)
l4<-layout.circle(network4)
l4<-layout.norm(l4, ymin=-1, ymax=1, xmin=-1, xmax=1)
l5<-layout.circle(network5)
l5<-layout.norm(l5, ymin=-1, ymax=1, xmin=-1, xmax=1)


pdf("figures/FigureXX.pdf")
layout(matrix(1:20, ncol=4, nrow=5))
par(oma=c(1, 1, 1, 1), mar=c(1, 1, 1, 1))

plot(1, xlim = c(min(dat1$Anjanaharidesud_low), max(dat1$Anjanaharidesud_high)),
  ylim = c(1, S1), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat1<-dat1[order(dat1$Anjanaharidesud_low),]
for (j in 1:S1) {
  lines(dat1[j,], c(j, j))
}

plot(1, xlim = c(min(dat2$Marojejy_low), max(dat2$Marojejy_high)),
     ylim = c(1, S2), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat2<-dat2[order(dat2$Marojejy_low),]
for (j in 1:S2) {
  lines(dat2[j,], c(j, j))
}

plot(1, xlim = c(min(dat3$Zahamena_low), max(dat3$Zahamena_high)),
     ylim = c(1, S3), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat3<-dat3[order(dat3$Zahamena_low),]
for (j in 1:S3) {
  lines(dat3[j,], c(j, j))
}

plot(1, xlim = c(min(dat4$Andringitra_low), max(dat4$Andringitra_high)),
     ylim = c(1, S4), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat4<-dat4[order(dat4$Andringitra_low),]
for (j in 1:S4) {
  lines(dat4[j,], c(j, j))
}

plot(1, xlim = c(min(dat5$Andohahela_low), max(dat5$Andohahela_high)),
     ylim = c(1, S5), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat5<-dat5[order(dat5$Andohahela_low),]
for (j in 1:S5) {
  lines(dat5[j,], c(j, j))
}

hist(log(range.sizes1), col="gray", border="gray", xlim=c(3,8), xlab="", 
     main="log(range size)", breaks=seq(from=3, to=8, by=0.3))

hist(log(range.sizes2), col="gray", border="gray", xlim=c(3,8), xlab="",
     main="", breaks=seq(from=3, to=8, by=0.3))

hist(log(range.sizes3), col="gray", border="gray", xlim=c(3,8), xlab="",
     main="", breaks=seq(from=3, to=8, by=0.3))

hist(log(range.sizes4), col="gray", border="gray", xlim=c(3,8), xlab="",
     main="", breaks=seq(from=3, to=8, by=0.3))

hist(log(range.sizes5), col="gray", border="gray", xlim=c(3,8), xlab="",
     main="", breaks=seq(from=3, to=8, by=0.3))

hist(log(overlap1), col="gray", border="gray", xlim=c(1,8), xlab="",
     main="log(range overlap)", breaks=seq(from=1, to=8, by=0.3))

hist(log(overlap2), col="gray", border="gray", xlim=c(1,8), xlab="",
     main="", breaks=seq(from=1, to=8, by=0.3))

hist(log(overlap3), col="gray", border="gray", xlim=c(1,8), xlab="",
     main="", breaks=seq(from=1, to=8, by=0.3))

hist(log(overlap4), col="gray", border="gray", xlim=c(1,8), xlab="",
     main="", breaks=seq(from=1, to=8, by=0.3))

hist(log(overlap5), col="gray", border="gray", xlim=c(1,8), xlab="",
     main="", breaks=seq(from=1, to=8, by=0.3))


plot(network1, rescale=F, layout=l1*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color=magma(5)[3], vertex.frame.color=magma(5)[3])

plot(network2, rescale=F, layout=l2*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color=magma(5)[3], vertex.frame.color=magma(5)[3])

plot(network3, rescale=F, layout=l3*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color=magma(5)[3], vertex.frame.color=magma(5)[3])

plot(network4, rescale=F, layout=l4*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color=magma(5)[3], vertex.frame.color=magma(5)[3])

plot(network5, rescale=F, layout=l5*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color=magma(5)[3], vertex.frame.color=magma(5)[3])

dev.off()
