rm(list=ls())

setwd("C:/Users/ferna/Dropbox/overnets")

source("source/functions.R")

x<-numeric()
for(i in 1:100) {
  xx<-rlnorm(n=1000)
  x[i]<-mean(xx)
}

S<-50

mat<-list()
res<-list()
for(i in 1:length(x)){
  mat[[i]]<-simRangeOver(x[i], S, TRUE)
  res[[i]]<-metRangeOver(mat[[i]])
}
