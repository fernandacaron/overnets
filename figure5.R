rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

library(sna)
library(igraph)

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

rho.MAH<-mean(log(range.sizes.MAH))/(max(dat_MAH$MAX)-min(dat_MAH$MIN))
rho.CDH<-mean(log(range.sizes.CDH))/(max(dat_CDH$MAX)-min(dat_CDH$MIN))
rho.CAN<-mean(log(range.sizes.CAN))/(max(dat_CAN$MAX)-min(dat_CAN$MIN))
rho.NAN<-mean(log(range.sizes.NAN))/(max(dat_NAN$MAX)-min(dat_NAN$MIN))
rho.TEP<-mean(log(range.sizes.TEP))/(max(dat_TEP$MAX)-min(dat_TEP$MIN))

diag(overlap.MAH)<-diag(overlap.CDH)<-diag(overlap.CAN)<-diag(overlap.NAN)<-
  diag(overlap.TEP)<-NA

over.MAH<-overlap.MAH[overlap.MAH > 0]
over.CDH<-overlap.CDH[overlap.CDH > 0]
over.CAN<-overlap.CAN[overlap.CAN > 0]
over.NAN<-overlap.NAN[overlap.NAN > 0]
over.TEP<-overlap.TEP[overlap.TEP > 0]

over.MAH<-na.omit(over.MAH)
over.CDH<-na.omit(over.CDH)
over.CAN<-na.omit(over.CAN)
over.NAN<-na.omit(over.NAN)
over.TEP<-na.omit(over.TEP)

net.MAH<-graph_from_adjacency_matrix(overlap.MAH, mode="undirected", weighted = T,
                                     diag=F)
net.CDH<-graph_from_adjacency_matrix(overlap.CDH, mode="undirected", weighted = T,
                                     diag=F)
net.CAN<-graph_from_adjacency_matrix(overlap.CAN, mode="undirected", weighted = T,
                                     diag=F)
net.NAN<-graph_from_adjacency_matrix(overlap.NAN, mode="undirected", weighted = T,
                                     diag=F)
net.TEP<-graph_from_adjacency_matrix(overlap.TEP, mode="undirected", weighted = T,
                                     diag=F)

l.MAH<-layout.circle(net.MAH)
l.MAH<-layout.norm(l.MAH, ymin=-1, ymax=1, xmin=-1, xmax=1)
l.CDH<-layout.circle(net.CDH)
l.CDH<-layout.norm(l.CDH, ymin=-1, ymax=1, xmin=-1, xmax=1)
l.CAN<-layout.circle(net.CAN)
l.CAN<-layout.norm(l.CAN, ymin=-1, ymax=1, xmin=-1, xmax=1)
l.NAN<-layout.circle(net.NAN)
l.NAN<-layout.norm(l.NAN, ymin=-1, ymax=1, xmin=-1, xmax=1)
l.TEP<-layout.circle(net.TEP)
l.TEP<-layout.norm(l.TEP, ymin=-1, ymax=1, xmin=-1, xmax=1)


pdf("figures/Figure5.1.pdf")
layout(matrix(1:20, ncol=4, nrow=5))
par(mar=c(4, 4, 1, 1))

plot(1, xlim = c(min(dat_MAH$MIN), max(dat_MAH$MAX)),
     ylim = c(1, S_MAH), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat_MAH<-dat_MAH[order(dat_MAH$MIN),]
for (j in 1:S_MAH) {
  lines(dat_MAH[j,], c(j, j), col=rgb(0,0,0,0.5))
}

plot(1, xlim = c(min(dat_CDH$MIN), max(dat_CDH$MAX)),
     ylim = c(1, S_CDH), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat_CDH<-dat_CDH[order(dat_CDH$MIN),]
for (j in 1:S_CDH) {
  lines(dat_CDH[j,], c(j, j), col=rgb(0,0,0,0.5))
}

plot(1, xlim = c(min(dat_CAN$MIN), max(dat_CAN$MAX)),
     ylim = c(1, S_CAN), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat_CAN<-dat_CAN[order(dat_CAN$MIN),]
for (j in 1:S_CAN) {
  lines(dat_CAN[j,], c(j, j), col=rgb(0,0,0,0.5))
}

plot(1, xlim = c(min(dat_NAN$MIN), max(dat_NAN$MAX)),
     ylim = c(1, S_NAN), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat_NAN<-dat_NAN[order(dat_NAN$MIN),]
for (j in 1:S_NAN) {
  lines(dat_NAN[j,], c(j, j), col=rgb(0,0,0,0.5))
}

plot(1, xlim = c(min(dat_TEP$MIN), max(dat_TEP$MAX)),
     ylim = c(1, S_TEP), type = "n", xlab = "Elevation", ylab = "Species ID", main = "")
dat_TEP<-dat_TEP[order(dat_TEP$MIN),]
for (j in 1:S_TEP) {
  lines(dat_TEP[j,], c(j, j), col=rgb(0,0,0,0.5))
}

hist(log(range.sizes.MAH), col="gray", border="gray", xlim=c(5,10), xlab="log(range size)", 
     main="", breaks=seq(from=5, to=9, by=0.2))

hist(log(range.sizes.CDH), col="gray", border="gray", xlim=c(4,10), xlab="log(range size)",
     main="", breaks=seq(from=5, to=9, by=0.2))

hist(log(range.sizes.CAN), col="gray", border="gray", xlim=c(3,10), xlab="log(range size)",
     main="", breaks=seq(from=4, to=9, by=0.2))

hist(log(range.sizes.NAN), col="gray", border="gray", xlim=c(3,10), xlab="log(range size)",
     main="", breaks=seq(from=4, to=9, by=0.2))

hist(log(range.sizes.TEP), col="gray", border="gray", xlim=c(5,10), xlab="log(range size)",
     main="", breaks=seq(from=5, to=9, by=0.2))

hist(log(over.MAH), col="gray", border="gray", xlim=c(2,10), xlab="log(range overlap)",
     main="", breaks=seq(from=3, to=9, by=0.2))

hist(log(over.CDH), col="gray", border="gray", xlim=c(2,10), xlab="log(range overlap)",
     main="", breaks=seq(from=3, to=9, by=0.2))

hist(log(over.CAN), col="gray", border="gray", xlim=c(2,10), xlab="log(range overlap)",
     main="", breaks=seq(from=3, to=9, by=0.2))

hist(log(over.NAN), col="gray", border="gray", xlim=c(2,10), xlab="log(range overlap)",
     main="", breaks=seq(from=3, to=9, by=0.2))

hist(log(over.TEP), col="gray", border="gray", xlim=c(2,10), xlab="log(range overlap)",
     main="", breaks=seq(from=3, to=9, by=0.2))


plot(net.MAH, rescale=F, layout=l.MAH*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color="black", vertex.frame.color="black")

plot(net.CDH, rescale=F, layout=l.CDH*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color="black", vertex.frame.color="black")

plot(net.CAN, rescale=F, layout=l.CAN*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color="black", vertex.frame.color="black")

plot(net.NAN, rescale=F, layout=l.NAN*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color="black", vertex.frame.color="black")

plot(net.TEP, rescale=F, layout=l.TEP*1, vertex.size=5, vertex.label=NA, edge.color="gray",
     vertex.color="black", vertex.frame.color="black")

dev.off()
