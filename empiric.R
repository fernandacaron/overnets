

setwd("C:/Users/ferna/Documents/IC/overnets")

source("codes/functions.R")

dat<-read.csv("data/birds_Anjanaharidesud_adjacency_matrix.csv", row.names=1)
dat<-na.omit(dat)
dat[,1]<-as.numeric(dat[,1])
dat[,2]<-as.numeric(dat[,2])
S<-nrow(dat)
  
overlap<-matrix(nrow=S, ncol=S)
for (i in 1:S) {
  for (j in 1:S) {
    pair <- rbind(dat[i, ], dat[j, ])
    overlap[i, j] <- overcalc(pair[1, ], pair[2, ])
  }
}

diag(overlap)<-NA
net<-network(overlap, directed = FALSE)

ggnet2(net, mode="circle", node.size=4, edge.color = 
         rgb(120,120,120, 0.2, maxColorValue = 255))

plot(
  1,
  xlim = c(860, 2000),
  ylim = c(1, S),
  type = "n",
  xlab = "D",
  ylab = "Species ID",
  main = ""
  
)
for (j in 1:S) {
  lines(dat[j,], c(j, j))
}
