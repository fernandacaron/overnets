rm(list=ls())

library(sna)
library(network)
library(ggplot2)
library(GGally)
library(egg)
library(viridis)

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

meanlog<-c(-3.5, -1.2)

S <- 10
range.sizes1 <- rlnorm(S, meanlog=meanlog[1])
#range.sizes2 <- rlnorm(S, meanlog=meanlog[2])
rho1<-round(mean(range.sizes1), digits = 2)
#rho2<-round(mean(range.sizes2), digits = 2)
rho1
#rho2
range.sizes1[range.sizes1>1]<-1
range.sizes2[range.sizes2>1]<-1
range.limits1 <- matrix(ncol = 2, nrow = S)
range.limits2 <- matrix(ncol = 2, nrow = S)
midpoint1 <- runif(S, min = 0, max = 1)
midpoint2 <- runif(S, min = 0, max = 1)
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
range.limits1[range.limits1>1]<-1
range.limits2[range.limits2>1]<-1

overlap1 <- overlap2 <- overlap3 <- matrix(nrow = S, ncol = S)
for (i in 1:S) {
  for (j in 1:S) {
    pair1 <- rbind(range.limits1[i, ], range.limits1[j, ])
    pair2 <- rbind(range.limits2[i, ], range.limits2[j, ])
    overlap1[i, j] <- overcalc(pair1[1, ], pair1[2, ])
    overlap2[i, j] <- overcalc(pair2[1, ], pair2[2, ])
  }
}

net1 <- network(overlap1, directed = FALSE)
net2 <- network(overlap2, directed = FALSE)
net3 <- network(overlap3, directed = FALSE)

network1<-graph_from_adjacency_matrix(overlap1, mode="undirected", weighted = T,
                                      diag=F)
network2<-graph_from_adjacency_matrix(overlap2, mode="undirected", weighted = T,
                                      diag=F)

l1<-layout.circle(network1)
l1 <- layout.norm(l1, ymin=-1, ymax=1, xmin=-1, xmax=1)
l2<-layout.circle(network2)
l2 <- layout.norm(l2, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf("figures/Figure1.pdf")
layout(matrix(1:4, ncol=2))
plot(1, xlim = c(0, 1), ylim = c(1, S), type = "n", xlab = "D", ylab = "Species ID",
     main = bquote(rho == 0.10))
for (j in 1:S) {
  lines(range.limits1[j,], c(j, j))
}
plot(1, xlim = c(0, 1), ylim = c(1, S), type = "n", xlab = "D", ylab = "Species ID",
     main = bquote(rho == 0.75))
for (j in 1:S) {
  lines(range.limits2[j,], c(j, j))
}
plot(network1, rescale=F, layout=l1*1.4, vertex.size=15, vertex.label=NA, edge.color="gray",
     vertex.color=viridis(3)[3], vertex.frame.color=viridis(3)[3])
plot(network2, rescale=F, layout=l2*1.4, vertex.size=15, vertex.label=NA, edge.color="gray",
     vertex.color=viridis(3)[1], vertex.frame.color=viridis(3)[1])

dev.off()

ggnet2(net1, node.color = viridis(3)[3], edge.color = 
            rgb(120,120,120, 1, maxColorValue = 255), mode = "circle")
ggnet2(net2, node.color = viridis(3)[1], edge.color = 
             rgb(120,120,120, 1, maxColorValue = 255), mode = "circle")
ggnet2(net3, node.color = viridis(3)[1], edge.color = 
         rgb(120,120,120, 1, maxColorValue = 255), mode = "circle")
