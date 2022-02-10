rm(list=ls())

library(network)
library(igraph)

source("functions.R")

#This code will simulate range sizes with two extreme values of rho to show 
#how it influences the range size distribution in the Domain and how influences
#the network structure

meanlog<-c(-3, -0.5)

S<-10

#Simulating range sizes for two different values of meanlog, calculating the 
#value of rho, and simulate range limits with these values, truncating when 
#larger than D.
range.sizes1 <- rlnorm(S, meanlog=meanlog[1])
range.sizes2 <- rlnorm(S, meanlog=meanlog[2])
rho1<-(mean(log(range.sizes1)))/log(exp(1))
rho2<-(mean(log(range.sizes2)))/log(exp(1))
range.limits1 <- matrix(ncol = 2, nrow = S)
range.limits2 <- matrix(ncol = 2, nrow = S)
midpoint1 <- runif(S, min = 0, max = exp(1))
midpoint2 <- runif(S, min = 0, max = exp(1))
for (i in 1:S) {
  half.rangesizes1 <- range.sizes1 / 2
  half.rangesizes2 <- range.sizes2 / 2
  range.limits1[i, 1] <- midpoint1[i] - half.rangesizes1[i]
  range.limits2[i, 1] <- midpoint2[i] - half.rangesizes2[i]
  range.limits1[i, 2] <- midpoint1[i] + half.rangesizes1[i]
  range.limits2[i, 2] <- midpoint2[i] + half.rangesizes2[i]
}
range.limits1[range.limits1<0]<-0
range.limits2[range.limits2<0]<-0
range.limits1[range.limits1>exp(1)]<-exp(1)
range.limits2[range.limits2>exp(1)]<-exp(1)

#Calculating the range overlap between species
range.limit.max1<-apply(range.limits1, 1, max)
range.limit.max2<-apply(range.limits2, 1, max)
range.limit.min1<-apply(range.limits1, 1, min)
range.limit.min2<-apply(range.limits2, 1, min)
overlap1<-overlap2<-matrix(nrow = S, ncol = S)
for(i in 1:S){
  pair.max1<-mapply(min, range.limit.max1[i], range.limit.max1)
  pair.max2<-mapply(min, range.limit.max2[i], range.limit.max2)
  pair.min1<-mapply(max, range.limit.min1[i], range.limit.min1)
  pair.min2<-mapply(max, range.limit.min2[i], range.limit.min2)
  overlap1[i,]<-pair.max1-pair.min1
  overlap2[i,]<-pair.max2-pair.min2
}
overlap1[overlap1 < 0] <- 0
overlap2[overlap2 < 0] <- 0
  
#Creating the network for each range overlap network
net1<-network(overlap1, directed = FALSE)
net2<-network(overlap2, directed = FALSE)

network1<-graph_from_adjacency_matrix(overlap1, mode="undirected", weighted = T,
                                      diag=F)
network2<-graph_from_adjacency_matrix(overlap2, mode="undirected", weighted = T,
                                      diag=F)

l1<-layout.circle(network1)
l1<-layout.norm(l1, ymin=-1, ymax=1, xmin=-1, xmax=1)
l2<-layout.circle(network2)
l2<-layout.norm(l2, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf("Figure2.pdf")
layout(matrix(1:4, ncol=2))
plot(1, xlim=c(0, exp(1)), ylim=c(1, S), type="n", xlab="D", ylab="Species ID",
     main=bquote(rho == -2.5))
for (j in 1:S) {
  lines(range.limits1[j,], c(j, j))
}
plot(1, xlim=c(0, exp(1)), ylim=c(1, S), type="n", xlab="D", ylab="Species ID",
     main=bquote(rho == -0.4))
for (j in 1:S) {
  lines(range.limits2[j,], c(j, j))
}
plot(network1, rescale=F, layout=l1*1.4, vertex.size=15, vertex.label=NA, 
     edge.color="gray", vertex.color="black", vertex.frame.color="black")
plot(network2, rescale=F, layout=l2*1.4, vertex.size=15, vertex.label=NA, 
     edge.color="gray", vertex.color="black", vertex.frame.color="black")

dev.off()
