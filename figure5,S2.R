rm(list=ls())

setwd("C:/Users/ferna/Documents/IC/overnets")

library(MASS)
library(sna)

source("codes/functions.R")

#This code will calculate range size for a empirical dataset, estimate the
#parameters of a fitted lognormal distribution to the range sizes, calculate
#the network metrics and simulate range overlap matrices to calculate the 
#expected values

dat<-read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv")

datMAH<-dat[dat$MAH=="Y",c(10,12)]
rownames(datMAH)<-paste(dat[dat$MAH=="Y",3], dat[dat$MAH=="Y",4], sep="_")
datMAH[,c(1:2)]<-apply(datMAH[,c(1:2)], 2, function(x) as.numeric((x)))

datCDH<-dat[dat$CDH=="Y",c(10,12)]
rownames(datCDH)<-paste(dat[dat$CDH=="Y",3], dat[dat$CDH=="Y",4], sep="_")
datCDH[,c(1:2)]<-apply(datCDH[,c(1:2)], 2, function(x) as.numeric((x)))

datCAN<-dat[dat$CAN=="Y",c(10,12)]
rownames(datCAN)<-paste(dat[dat$CAN=="Y",3], dat[dat$CAN=="Y",4], sep="_")
datCAN[,c(1:2)]<-apply(datCAN[,c(1:2)], 2, function(x) as.numeric((x)))

datNAN<-dat[dat$NAN=="Y",c(10,12)]
rownames(datNAN)<-paste(dat[dat$NAN=="Y",3], dat[dat$NAN=="Y",4], sep="_")
datNAN[,c(1:2)]<-apply(datNAN[,c(1:2)], 2, function(x) as.numeric((x)))

datTEP<-dat[dat$TEP=="Y",c(10,12)]
rownames(datTEP)<-paste(dat[dat$TEP=="Y",3], dat[dat$TEP=="Y",4], sep="_")
datTEP[,c(1:2)]<-apply(datTEP[,c(1:2)], 2, function(x) as.numeric((x)))

datMAH<-na.omit(datMAH)
datCDH<-na.omit(datCDH)
datCAN<-na.omit(datCAN)
datNAN<-na.omit(datNAN)
datTEP<-na.omit(datTEP)

range.sizes.MAH<-numeric()
for(i in 1:nrow(datMAH)){
  range.sizes.MAH[i]<-datMAH$MAX[i]-datMAH$MIN[i]
}
names(range.sizes.MAH)<-rownames(datMAH)
range.sizes.MAH<-range.sizes.MAH[range.sizes.MAH>0]

range.sizes.CDH<-numeric()
for(i in 1:nrow(datCDH)){
  range.sizes.CDH[i]<-datCDH$MAX[i]-datCDH$MIN[i]
}
names(range.sizes.CDH)<-rownames(datCDH)
range.sizes.CDH<-range.sizes.CDH[range.sizes.CDH>0]

range.sizes.CAN<-numeric()
for(i in 1:nrow(datCAN)){
  range.sizes.CAN[i]<-datCAN$MAX[i]-datCAN$MIN[i]
}
names(range.sizes.CAN)<-rownames(datCAN)
range.sizes.CAN<-range.sizes.CAN[range.sizes.CAN>0]

range.sizes.NAN<-numeric()
for(i in 1:nrow(datNAN)){
  range.sizes.NAN[i]<-datNAN$MAX[i]-datNAN$MIN[i]
}
names(range.sizes.NAN)<-rownames(datNAN)
range.sizes.NAN<-range.sizes.NAN[range.sizes.NAN>0]

range.sizes.TEP<-numeric()
for(i in 1:nrow(datTEP)){
  range.sizes.TEP[i]<-datTEP$MAX[i]-datTEP$MIN[i]
}
names(range.sizes.TEP)<-rownames(datTEP)
range.sizes.TEP<-range.sizes.TEP[range.sizes.TEP>0]

datMAH<-datMAH[rownames(datMAH) %in% names(range.sizes.MAH),]
datCDH<-datCDH[rownames(datCDH) %in% names(range.sizes.CDH),]
datCAN<-datCAN[rownames(datCAN) %in% names(range.sizes.CAN),]
datNAN<-datNAN[rownames(datNAN) %in% names(range.sizes.NAN),]
datTEP<-datTEP[rownames(datTEP) %in% names(range.sizes.TEP),]

S.MAH<-nrow(datMAH)
S.CDH<-nrow(datCDH)
S.CAN<-nrow(datCAN)
S.NAN<-nrow(datNAN)
S.TEP<-nrow(datTEP)

overlap.MAH<-overlap.matrix(datMAH)
overlap.CDH<-overlap.matrix(datCDH)
overlap.CAN<-overlap.matrix(datCAN)
overlap.NAN<-overlap.matrix(datNAN)
overlap.TEP<-overlap.matrix(datTEP)

rho.MAH<-mean(log(range.sizes.MAH))/log((max(datMAH$MAX)-min(datMAH$MIN)))
rho.CDH<-mean(log(range.sizes.CDH))/log((max(datCDH$MAX)-min(datCDH$MIN)))
rho.CAN<-mean(log(range.sizes.CAN))/log((max(datCAN$MAX)-min(datCAN$MIN)))
rho.NAN<-mean(log(range.sizes.NAN))/log((max(datNAN$MAX)-min(datNAN$MIN)))
rho.TEP<-mean(log(range.sizes.TEP))/log((max(datTEP$MAX)-min(datTEP$MIN)))

over.MAH<-overlap.MAH
over.CDH<-overlap.CDH
over.CAN<-overlap.CAN
over.NAN<-overlap.NAN
over.TEP<-overlap.TEP

diag(over.MAH)<-diag(over.CDH)<-diag(over.CAN)<-diag(over.NAN)<-
  diag(over.TEP)<-0

over.MAH[over.MAH>0]<-1
over.CDH[over.CDH>0]<-1
over.CAN[over.CAN>0]<-1
over.NAN[over.NAN>0]<-1
over.TEP[over.TEP>0]<-1

overlap.MAH<-overlap.MAH[overlap.MAH>0]
overlap.CDH<-overlap.CDH[overlap.CDH>0]
overlap.CAN<-overlap.CAN[overlap.CAN>0]
overlap.NAN<-overlap.NAN[overlap.NAN>0]
overlap.TEP<-overlap.TEP[overlap.TEP>0]

meanNumOverSpp.MAH<-mean(rowSums(over.MAH))
meanNumOverSpp.CDH<-mean(rowSums(over.CDH))
meanNumOverSpp.CAN<-mean(rowSums(over.CAN))
meanNumOverSpp.NAN<-mean(rowSums(over.NAN))
meanNumOverSpp.TEP<-mean(rowSums(over.TEP))

meanlog.MAH<-fitdistr(range.sizes.MAH,"lognormal")$estimate[[1]]
meanlog.CDH<-fitdistr(range.sizes.CDH,"lognormal")$estimate[[1]]
meanlog.CAN<-fitdistr(range.sizes.CAN,"lognormal")$estimate[[1]]
meanlog.NAN<-fitdistr(range.sizes.NAN,"lognormal")$estimate[[1]]
meanlog.TEP<-fitdistr(range.sizes.TEP,"lognormal")$estimate[[1]]

sdlog.MAH<-fitdistr(range.sizes.MAH,"lognormal")$estimate[[2]]
sdlog.CDH<-fitdistr(range.sizes.CDH,"lognormal")$estimate[[2]]
sdlog.CAN<-fitdistr(range.sizes.CAN,"lognormal")$estimate[[2]]
sdlog.NAN<-fitdistr(range.sizes.NAN,"lognormal")$estimate[[2]]
sdlog.TEP<-fitdistr(range.sizes.TEP,"lognormal")$estimate[[2]]

netMAH<-matrix(nrow=S.MAH, ncol=3)
netCDH<-matrix(nrow=S.CDH, ncol=3)
netCAN<-matrix(nrow=S.CAN, ncol=3)
netNAN<-matrix(nrow=S.NAN, ncol=3)
netTEP<-matrix(nrow=S.TEP, ncol=3)

netMAH[,1]<-sort(degree(over.MAH, gmode="graph"), decreasing=T)
netMAH[,2]<-sort(betweenness(over.MAH, gmode="graph"), decreasing=T)
netMAH[,3]<-sort(closeness(over.MAH, gmode="graph"), decreasing=T)

netCDH[,1]<-sort(degree(over.CDH, gmode="graph"), decreasing=T)
netCDH[,2]<-sort(betweenness(over.CDH, gmode="graph"), decreasing=T)
netCDH[,3]<-sort(closeness(over.CDH, gmode="graph"), decreasing=T)

netCAN[,1]<-sort(degree(over.CAN, gmode="graph"), decreasing=T)
netCAN[,2]<-sort(betweenness(over.CAN, gmode="graph"), decreasing=T)
netCAN[,3]<-sort(closeness(over.CAN, gmode="graph"), decreasing=T)

netNAN[,1]<-sort(degree(over.NAN, gmode="graph"), decreasing=T)
netNAN[,2]<-sort(betweenness(over.NAN, gmode="graph"), decreasing=T)
netNAN[,3]<-sort(closeness(over.NAN, gmode="graph"), decreasing=T)

netTEP[,1]<-sort(degree(over.TEP, gmode="graph"), decreasing=T)
netTEP[,2]<-sort(betweenness(over.TEP, gmode="graph"), decreasing=T)
netTEP[,3]<-sort(closeness(over.TEP, gmode="graph"), decreasing=T)

nsim<-100
simMeanNumOverSpp.MAH.NT<-simRho.MAH.NT<-simMeanOver.MAH.NT<-simVarOver.MAH.NT<-
  numeric()
simOver.MAH.NT<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.MAH, sdlog=sdlog.MAH, 
                                       S=S.MAH, truncate=F, Dmin=min(datMAH$MIN), 
                                       Dmax=max(datMAH$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.MAH.NT[i]<-mean(rowSums(overlap))
  simRho.MAH.NT[i]<-dat$rho
  
  simOver.MAH.NT[[i]]<-overlap

  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.MAH.NT[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.MAH.NT[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.MAH.NT<-simBet.MAH.NT<-simClo.MAH.NT<-matrix(nrow=S.MAH, ncol=nsim)
for(i in 1:nsim){
  simDeg.MAH.NT[,i]<-sort(degree(simOver.MAH.NT[[i]], gmode="graph"), decreasing=T)
  simBet.MAH.NT[,i]<-sort(betweenness(simOver.MAH.NT[[i]], gmode="graph"), decreasing=T)
  simClo.MAH.NT[,i]<-sort(closeness(simOver.MAH.NT[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.CDH.NT<-simRho.CDH.NT<-simMeanOver.CDH.NT<-simVarOver.CDH.NT<-
  numeric()
simOver.CDH.NT<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.CDH, sdlog=sdlog.CDH, 
                    S=S.CDH, truncate=F, Dmin=min(datCDH$MIN), 
                    Dmax=max(datCDH$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.CDH.NT[i]<-mean(rowSums(overlap))
  simRho.CDH.NT[i]<-dat$rho
  
  simOver.CDH.NT[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.CDH.NT[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.CDH.NT[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.CDH.NT<-simBet.CDH.NT<-simClo.CDH.NT<-matrix(nrow=S.CDH, ncol=nsim)
for(i in 1:nsim){
  simDeg.CDH.NT[,i]<-sort(degree(simOver.CDH.NT[[i]], gmode="graph"), decreasing=T)
  simBet.CDH.NT[,i]<-sort(betweenness(simOver.CDH.NT[[i]], gmode="graph"), decreasing=T)
  simClo.CDH.NT[,i]<-sort(closeness(simOver.CDH.NT[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.CAN.NT<-simRho.CAN.NT<-simMeanOver.CAN.NT<-simVarOver.CAN.NT<-
  numeric()
simOver.CAN.NT<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.CAN, sdlog=sdlog.CAN, 
                    S=S.CAN, truncate=F, Dmin=min(datCAN$MIN), 
                    Dmax=max(datCAN$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.CAN.NT[i]<-mean(rowSums(overlap))
  simRho.CAN.NT[i]<-dat$rho
  
  simOver.CAN.NT[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.CAN.NT[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.CAN.NT[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.CAN.NT<-simBet.CAN.NT<-simClo.CAN.NT<-matrix(nrow=S.CAN, ncol=nsim)
for(i in 1:nsim){
  simDeg.CAN.NT[,i]<-sort(degree(simOver.CAN.NT[[i]], gmode="graph"), decreasing=T)
  simBet.CAN.NT[,i]<-sort(betweenness(simOver.CAN.NT[[i]], gmode="graph"), decreasing=T)
  simClo.CAN.NT[,i]<-sort(closeness(simOver.CAN.NT[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.NAN.NT<-simRho.NAN.NT<-simMeanOver.NAN.NT<-simVarOver.NAN.NT<-
  numeric()
simOver.NAN.NT<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.NAN, sdlog=sdlog.NAN, 
                    S=S.NAN, truncate=F, Dmin=min(datNAN$MIN), 
                    Dmax=max(datNAN$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.NAN.NT[i]<-mean(rowSums(overlap))
  simRho.NAN.NT[i]<-dat$rho
  
  simOver.NAN.NT[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.NAN.NT[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.NAN.NT[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.NAN.NT<-simBet.NAN.NT<-simClo.NAN.NT<-matrix(nrow=S.NAN, ncol=nsim)
for(i in 1:nsim){
  simDeg.NAN.NT[,i]<-sort(degree(simOver.NAN.NT[[i]], gmode="graph"), decreasing=T)
  simBet.NAN.NT[,i]<-sort(betweenness(simOver.NAN.NT[[i]], gmode="graph"), decreasing=T)
  simClo.NAN.NT[,i]<-sort(closeness(simOver.NAN.NT[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.TEP.NT<-simRho.TEP.NT<-simMeanOver.TEP.NT<-simVarOver.TEP.NT<-
  numeric()
simOver.TEP.NT<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.TEP, sdlog=sdlog.TEP, 
                    S=S.TEP, truncate=F, Dmin=min(datTEP$MIN), 
                    Dmax=max(datTEP$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.TEP.NT[i]<-mean(rowSums(overlap))
  simRho.TEP.NT[i]<-dat$rho
  
  simOver.TEP.NT[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.TEP.NT[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.TEP.NT[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.TEP.NT<-simBet.TEP.NT<-simClo.TEP.NT<-matrix(nrow=S.TEP, ncol=nsim)
for(i in 1:nsim){
  simDeg.TEP.NT[,i]<-sort(degree(simOver.TEP.NT[[i]], gmode="graph"), decreasing=T)
  simBet.TEP.NT[,i]<-sort(betweenness(simOver.TEP.NT[[i]], gmode="graph"), decreasing=T)
  simClo.TEP.NT[,i]<-sort(closeness(simOver.TEP.NT[[i]], gmode="graph"), decreasing=T)
}

pdf("figures/Figure5.pdf", width=12, height=14)
layout(matrix(1:35, nrow=7, ncol=5, byrow=T))
par(mar=c(4, 4, 1, 1))

col1<-rgb(33/255,12/255,74/255, 0.8)
col2<-rgb(137/255,34/255,106/255, 0.8)
col3<-rgb(187/255,55/255,84/255, 0.8)
col4<-rgb(227/255,89/255,50/255, 0.8)
col5<-rgb(249/255,140/255,10/255, 0.8)

#Linha 1: gráfico com linhas ordenadas
plot(1, xlim = c(min(datMAH$MIN), max(datMAH$MAX)),
     ylim = c(1, S.MAH), type = "n", xlab = "", ylab = "Species ID", main = "")
datMAH<-datMAH[order(datMAH$MIN),]
for (j in 1:S.MAH) {
  lines(datMAH[j,], c(j, j), col=col1)
}

plot(1, xlim = c(min(datCDH$MIN), max(datCDH$MAX)),
     ylim = c(1, S.CDH), type = "n", xlab = "", ylab = "", main = "")
datCDH<-datCDH[order(datCDH$MIN),]
for (j in 1:S.CDH) {
  lines(datCDH[j,], c(j, j), col=col2)
}

plot(1, xlim = c(min(datCAN$MIN), max(datCAN$MAX)),
     ylim = c(1, S.CAN), type = "n", xlab = "Elevation", ylab = "", main = "")
datCAN<-datCAN[order(datCAN$MIN),]
for (j in 1:S.CAN) {
  lines(datCAN[j,], c(j, j), col=col3)
}

plot(1, xlim = c(min(datNAN$MIN), max(datNAN$MAX)),
     ylim = c(1, S.NAN), type = "n", xlab = "", ylab = "", main = "")
datNAN<-datNAN[order(datNAN$MIN),]
for (j in 1:S.NAN) {
  lines(datNAN[j,], c(j, j), col=col4)
}

plot(1, xlim = c(min(datTEP$MIN), max(datTEP$MAX)),
     ylim = c(1, S.TEP), type = "n", xlab = "", ylab = "", main = "")
datTEP<-datTEP[order(datTEP$MIN),]
for (j in 1:S.TEP) {
  lines(datTEP[j,], c(j, j), col=col5)
}

#Linha 2: Range size distribution
hist(log(range.sizes.MAH), col=col1, border=col1, 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

hist(log(range.sizes.CDH), col=col2, border=col2, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

hist(log(range.sizes.CAN), col=col3, border=col3, ylab="",  
     xlim=c(2,10), xlab="Range size", main="", breaks=seq(from=4, to=9, by=0.15))

hist(log(range.sizes.NAN), col=col4, border=col4, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=4, to=9, by=0.15))

hist(log(range.sizes.TEP), col=col5, border=col5, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

#Linha 3: Range overlap distribution e a média de várias reps dos simulados

hist(log(overlap.MAH), col=col1, border=col1, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), freq=F, xlim=c(2,10))
hist(simMeanOver.MAH.NT, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T, 
     border=rgb(0.7,0.7,0.7,0.6),
     breaks=seq(from=min(simMeanOver.MAH.NT)-0.1, 
                to=max(simMeanOver.MAH.NT)+0.3, by=0.15))

hist(log(overlap.CDH), col=col2, border=col2, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.CDH.NT, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T,
     border=rgb(0.7,0.7,0.7,0.6), 
     breaks=seq(from=min(simMeanOver.CDH.NT)-0.1, 
                to=max(simMeanOver.CDH.NT)+0.3, by=0.15))

hist(log(overlap.CAN), col=col3, border=col3, main="", xlab="Range overlap",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.CAN.NT, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T,
     border=rgb(0.7,0.7,0.7,0.6), 
     breaks=seq(from=min(simMeanOver.CAN.NT)-0.1, 
                to=max(simMeanOver.CAN.NT)+0.3, by=0.15))

hist(log(overlap.NAN), col=col4, border=col4, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.NAN.NT, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.6),
     breaks=seq(from=min(simMeanOver.NAN.NT)-0.1, 
                to=max(simMeanOver.NAN.NT)+0.3, by=0.15))

hist(log(overlap.TEP), col=col5, border=col5, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.TEP.NT, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T, 
     border=rgb(0.7,0.7,0.7,0.6),  
     breaks=seq(from=min(simMeanOver.TEP.NT)-0.1, 
                to=max(simMeanOver.TEP.NT)+0.3, by=0.15))

#Linha 4: Distribuição de mean number of overlapping species 
###(linha do observado + distribuição dos simulados)
hist(simMeanNumOverSpp.MAH.NT, col=col1, border=col1, 
     xlab="", main="", xlim=c(200,400),
     breaks=seq(from=min(simMeanNumOverSpp.MAH.NT)-10, 
                to=max(simMeanNumOverSpp.MAH.NT)+10, by=5))
abline(v=meanNumOverSpp.MAH, col=col1)

hist(simMeanNumOverSpp.CDH.NT, col=col2, border=col2, ylab="",
     xlab="", main="", xlim=c(140,270),
     breaks=seq(from=min(simMeanNumOverSpp.CDH.NT)-10, 
                to=max(simMeanNumOverSpp.CDH.NT)+10, by=3.5))
abline(v=meanNumOverSpp.CDH, col=col2)

hist(simMeanNumOverSpp.CAN.NT, col=col3, border=col3, ylab="",
     xlab="Average number of overlapping species", main="",
     xlim=c(350,600),
     breaks=seq(from=min(simMeanNumOverSpp.CAN.NT)-10, 
                to=max(simMeanNumOverSpp.CAN.NT)+10, by=6.3))
abline(v=meanNumOverSpp.CAN, col=col3)

hist(simMeanNumOverSpp.NAN.NT, col=col4, border=col4, ylab="",
     xlab="", main="", xlim=c(350,600),
     breaks=seq(from=min(simMeanNumOverSpp.NAN.NT)-10, 
                to=max(simMeanNumOverSpp.NAN.NT)+10, by=6.3))
abline(v=meanNumOverSpp.NAN, col=col4)

hist(simMeanNumOverSpp.TEP.NT, col=col5, border=col5, ylab="",
     xlab="", main="", xlim=c(60,140),
     breaks=seq(from=min(simMeanNumOverSpp.TEP.NT)-10, 
                to=max(simMeanNumOverSpp.TEP.NT)+10, by=2))
abline(v=meanNumOverSpp.TEP, col=col5)


#Linha 5: Degree distribution

plot(simDeg.MAH.NT[,1], main="", ylab="Degree", xlab="", type="l", 
     ylim=range(simDeg.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,1], col=col1) 

plot(simDeg.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,1], col=col2) 

plot(simDeg.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simDeg.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,1], col=col3) 

plot(simDeg.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,1], col=col4) 

plot(simDeg.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,1], col=col5)

#Linha 6: Betweenness distribution.

plot(simBet.MAH.NT[,1], main="", ylab="Betweenness", xlab="", type="l", 
     ylim=range(simBet.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,2], col=col1) 

plot(simBet.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,2], col=col2) 

plot(simBet.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simBet.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,2], col=col3) 

plot(simBet.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,2], col=col4) 

plot(simBet.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,2], col=col5)

#Linha 7: Closeness distribution.

plot(simClo.MAH.NT[,1], main="", ylab="Closeness", xlab="", type="l", 
     ylim=range(simClo.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,3], col=col1) 

plot(simClo.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,3], col=col2) 

plot(simClo.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simClo.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,3], col=col3) 

plot(simClo.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,3], col=col4) 

plot(simClo.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,3], col=col5)

dev.off()

nsim<-100
simMeanNumOverSpp.MAH.T<-simRho.MAH.T<-simMeanOver.MAH.T<-simVarOver.MAH.T<-
  numeric()
simOver.MAH.T<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.MAH, sdlog=sdlog.MAH, 
                    S=S.MAH, truncate=T, Dmin=min(datMAH$MIN), 
                    Dmax=max(datMAH$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.MAH.T[i]<-mean(rowSums(overlap))
  simRho.MAH.T[i]<-dat$rho
  
  simOver.MAH.T[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.MAH.T[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.MAH.T[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.MAH.T<-simBet.MAH.T<-simClo.MAH.T<-matrix(nrow=S.MAH, ncol=nsim)
for(i in 1:nsim){
  simDeg.MAH.T[,i]<-sort(degree(simOver.MAH.T[[i]], gmode="graph"), decreasing=T)
  simBet.MAH.T[,i]<-sort(betweenness(simOver.MAH.T[[i]], gmode="graph"), decreasing=T)
  simClo.MAH.T[,i]<-sort(closeness(simOver.MAH.T[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.CDH.T<-simRho.CDH.T<-simMeanOver.CDH.T<-simVarOver.CDH.T<-
  numeric()
simOver.CDH.T<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.CDH, sdlog=sdlog.CDH, 
                    S=S.CDH, truncate=T, Dmin=min(datCDH$MIN), 
                    Dmax=max(datCDH$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.CDH.T[i]<-mean(rowSums(overlap))
  simRho.CDH.T[i]<-dat$rho
  
  simOver.CDH.T[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.CDH.T[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.CDH.T[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.CDH.T<-simBet.CDH.T<-simClo.CDH.T<-matrix(nrow=S.CDH, ncol=nsim)
for(i in 1:nsim){
  simDeg.CDH.T[,i]<-sort(degree(simOver.CDH.T[[i]], gmode="graph"), decreasing=T)
  simBet.CDH.T[,i]<-sort(betweenness(simOver.CDH.T[[i]], gmode="graph"), decreasing=T)
  simClo.CDH.T[,i]<-sort(closeness(simOver.CDH.T[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.CAN.T<-simRho.CAN.T<-simMeanOver.CAN.T<-simVarOver.CAN.T<-
  numeric()
simOver.CAN.T<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.CAN, sdlog=sdlog.CAN, 
                    S=S.CAN, truncate=T, Dmin=min(datCAN$MIN), 
                    Dmax=max(datCAN$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.CAN.T[i]<-mean(rowSums(overlap))
  simRho.CAN.T[i]<-dat$rho
  
  simOver.CAN.T[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.CAN.T[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.CAN.T[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.CAN.T<-simBet.CAN.T<-simClo.CAN.T<-matrix(nrow=S.CAN, ncol=nsim)
for(i in 1:nsim){
  simDeg.CAN.T[,i]<-sort(degree(simOver.CAN.T[[i]], gmode="graph"), decreasing=T)
  simBet.CAN.T[,i]<-sort(betweenness(simOver.CAN.T[[i]], gmode="graph"), decreasing=T)
  simClo.CAN.T[,i]<-sort(closeness(simOver.CAN.T[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.NAN.T<-simRho.NAN.T<-simMeanOver.NAN.T<-simVarOver.NAN.T<-
  numeric()
simOver.NAN.T<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.NAN, sdlog=sdlog.NAN, 
                    S=S.NAN, truncate=T, Dmin=min(datNAN$MIN), 
                    Dmax=max(datNAN$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.NAN.T[i]<-mean(rowSums(overlap))
  simRho.NAN.T[i]<-dat$rho
  
  simOver.NAN.T[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.NAN.T[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.NAN.T[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.NAN.T<-simBet.NAN.T<-simClo.NAN.T<-matrix(nrow=S.NAN, ncol=nsim)
for(i in 1:nsim){
  simDeg.NAN.T[,i]<-sort(degree(simOver.NAN.T[[i]], gmode="graph"), decreasing=T)
  simBet.NAN.T[,i]<-sort(betweenness(simOver.NAN.T[[i]], gmode="graph"), decreasing=T)
  simClo.NAN.T[,i]<-sort(closeness(simOver.NAN.T[[i]], gmode="graph"), decreasing=T)
}

simMeanNumOverSpp.TEP.T<-simRho.TEP.T<-simMeanOver.TEP.T<-simVarOver.TEP.T<-
  numeric()
simOver.TEP.T<-list()
for(i in 1:1000){
  dat<-simRangeOver(meanlog=meanlog.TEP, sdlog=sdlog.TEP, 
                    S=S.TEP, truncate=T, Dmin=min(datTEP$MIN), 
                    Dmax=max(datTEP$MAX)) 
  overlap<-dat$overlap
  over<-dat$overlap
  diag(overlap)<-0
  overlap[overlap>0]<-1
  simMeanNumOverSpp.TEP.T[i]<-mean(rowSums(overlap))
  simRho.TEP.T[i]<-dat$rho
  
  simOver.TEP.T[[i]]<-overlap
  
  diag(over)<-NA
  over<-over[over>0]
  simMeanOver.TEP.T[i]<-mean(log(over), na.rm = TRUE)
  simVarOver.TEP.T[i]<-var(as.vector(log(over)), na.rm = TRUE)
}

simDeg.TEP.T<-simBet.TEP.T<-simClo.TEP.T<-matrix(nrow=S.TEP, ncol=nsim)
for(i in 1:nsim){
  simDeg.TEP.T[,i]<-sort(degree(simOver.TEP.T[[i]], gmode="graph"), decreasing=T)
  simBet.TEP.T[,i]<-sort(betweenness(simOver.TEP.T[[i]], gmode="graph"), decreasing=T)
  simClo.TEP.T[,i]<-sort(closeness(simOver.TEP.T[[i]], gmode="graph"), decreasing=T)
}

pdf("figures/FigureS2.pdf", width=12, height=14)
layout(matrix(1:35, nrow=7, ncol=5, byrow=T))
par(mar=c(4, 4, 1, 1))

col1<-rgb(33/255,12/255,74/255, 0.8)
col2<-rgb(137/255,34/255,106/255, 0.8)
col3<-rgb(187/255,55/255,84/255, 0.8)
col4<-rgb(227/255,89/255,50/255, 0.8)
col5<-rgb(249/255,140/255,10/255, 0.8)

#Linha 1: gráfico com linhas ordenadas
plot(1, xlim = c(min(datMAH$MIN), max(datMAH$MAX)),
     ylim = c(1, S.MAH), type = "n", xlab = "", ylab = "Species ID", main = "")
datMAH<-datMAH[order(datMAH$MIN),]
for (j in 1:S.MAH) {
  lines(datMAH[j,], c(j, j), col=col1)
}

plot(1, xlim = c(min(datCDH$MIN), max(datCDH$MAX)),
     ylim = c(1, S.CDH), type = "n", xlab = "", ylab = "", main = "")
datCDH<-datCDH[order(datCDH$MIN),]
for (j in 1:S.CDH) {
  lines(datCDH[j,], c(j, j), col=col2)
}

plot(1, xlim = c(min(datCAN$MIN), max(datCAN$MAX)),
     ylim = c(1, S.CAN), type = "n", xlab = "Elevation", ylab = "", main = "")
datCAN<-datCAN[order(datCAN$MIN),]
for (j in 1:S.CAN) {
  lines(datCAN[j,], c(j, j), col=col3)
}

plot(1, xlim = c(min(datNAN$MIN), max(datNAN$MAX)),
     ylim = c(1, S.NAN), type = "n", xlab = "", ylab = "", main = "")
datNAN<-datNAN[order(datNAN$MIN),]
for (j in 1:S.NAN) {
  lines(datNAN[j,], c(j, j), col=col4)
}

plot(1, xlim = c(min(datTEP$MIN), max(datTEP$MAX)),
     ylim = c(1, S.TEP), type = "n", xlab = "", ylab = "", main = "")
datTEP<-datTEP[order(datTEP$MIN),]
for (j in 1:S.TEP) {
  lines(datTEP[j,], c(j, j), col=col5)
}

#Linha 2: Range size distribution
hist(log(range.sizes.MAH), col=col1, border=col1, 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

hist(log(range.sizes.CDH), col=col2, border=col2, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

hist(log(range.sizes.CAN), col=col3, border=col3, ylab="",  
     xlim=c(2,10), xlab="Range size", main="", breaks=seq(from=4, to=9, by=0.15))

hist(log(range.sizes.NAN), col=col4, border=col4, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=4, to=9, by=0.15))

hist(log(range.sizes.TEP), col=col5, border=col5, ylab="", 
     xlim=c(2,10), xlab="", main="", breaks=seq(from=5, to=9, by=0.15))

#Linha 3: Range overlap distribution e a média de várias reps dos simulados

hist(log(overlap.MAH), col=col1, border=col1, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), freq=F, xlim=c(2,10))
hist(simMeanOver.MAH.T, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T, 
     border=rgb(0.7,0.7,0.7,0.6),
     breaks=seq(from=min(simMeanOver.MAH.T)-0.1, 
                to=max(simMeanOver.MAH.T)+0.3, by=0.15))

hist(log(overlap.CDH), col=col2, border=col2, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.CDH.T, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T,
     border=rgb(0.7,0.7,0.7,0.6), 
     breaks=seq(from=min(simMeanOver.CDH.T)-0.1, 
                to=max(simMeanOver.CDH.T)+0.3, by=0.15))

hist(log(overlap.CAN), col=col3, border=col3, main="", xlab="Range overlap",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.CAN.T, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T,
     border=rgb(0.7,0.7,0.7,0.6), 
     breaks=seq(from=min(simMeanOver.CAN.T)-0.1, 
                to=max(simMeanOver.CAN.T)+0.3, by=0.15))

hist(log(overlap.NAN), col=col4, border=col4, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.NAN.T, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.6),
     breaks=seq(from=min(simMeanOver.NAN.T)-0.1, 
                to=max(simMeanOver.NAN.T)+0.3, by=0.15))

hist(log(overlap.TEP), col=col5, border=col5, main="", xlab="",
     breaks=seq(from=3, to=9, by=0.15), ylab="", freq=F, xlim=c(2,10))
hist(simMeanOver.TEP.T, col=rgb(0.7,0.7,0.7,0.6), add=T, freq=T, 
     border=rgb(0.7,0.7,0.7,0.6),  
     breaks=seq(from=min(simMeanOver.TEP.T)-0.1, 
                to=max(simMeanOver.TEP.T)+0.3, by=0.15))

#Linha 4: Distribuição de mean number of overlapping species 
###(linha do observado + distribuição dos simulados)
hist(simMeanNumOverSpp.MAH.T, col=col1, border=col1, 
     xlab="", main="", xlim=c(200,400),
     breaks=seq(from=min(simMeanNumOverSpp.MAH.T)-10, 
                to=max(simMeanNumOverSpp.MAH.T)+10, by=5))
abline(v=meanNumOverSpp.MAH, col=col1)

hist(simMeanNumOverSpp.CDH.T, col=col2, border=col2, ylab="",
     xlab="", main="", xlim=c(140,270),
     breaks=seq(from=min(simMeanNumOverSpp.CDH.T)-10, 
                to=max(simMeanNumOverSpp.CDH.T)+10, by=3.5))
abline(v=meanNumOverSpp.CDH, col=col2)

hist(simMeanNumOverSpp.CAN.T, col=col3, border=col3, ylab="",
     xlab="Average number of overlapping species", main="",
     xlim=c(350,600),
     breaks=seq(from=min(simMeanNumOverSpp.CAN.T)-10, 
                to=max(simMeanNumOverSpp.CAN.T)+10, by=6.3))
abline(v=meanNumOverSpp.CAN, col=col3)

hist(simMeanNumOverSpp.NAN.T, col=col4, border=col4, ylab="",
     xlab="", main="", xlim=c(350,600),
     breaks=seq(from=min(simMeanNumOverSpp.NAN.T)-10, 
                to=max(simMeanNumOverSpp.NAN.T)+10, by=6.3))
abline(v=meanNumOverSpp.NAN, col=col4)

hist(simMeanNumOverSpp.TEP.T, col=col5, border=col5, ylab="",
     xlab="", main="", xlim=c(60,140),
     breaks=seq(from=min(simMeanNumOverSpp.TEP.T)-10, 
                to=max(simMeanNumOverSpp.TEP.T)+10, by=2))
abline(v=meanNumOverSpp.TEP, col=col5)


#Linha 5: Degree distribution

plot(simDeg.MAH.T[,1], main="", ylab="Degree", xlab="", type="l", 
     ylim=range(simDeg.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,1], col=col1) 

plot(simDeg.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,1], col=col2) 

plot(simDeg.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simDeg.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,1], col=col3) 

plot(simDeg.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,1], col=col4) 

plot(simDeg.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simDeg.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,1], col=col5)

#Linha 6: Betweenness distribution.

plot(simBet.MAH.T[,1], main="", ylab="Betweenness", xlab="", type="l", 
     ylim=range(simBet.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,2], col=col1) 

plot(simBet.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,2], col=col2) 

plot(simBet.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simBet.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,2], col=col3) 

plot(simBet.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,2], col=col4) 

plot(simBet.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simBet.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,2], col=col5)

#Linha 7: Closeness distribution.

plot(simClo.MAH.T[,1], main="", ylab="Closeness", xlab="", type="l", 
     ylim=range(simClo.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,3], col=col1) 

plot(simClo.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,3], col=col2) 

plot(simClo.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simClo.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,3], col=col3) 

plot(simClo.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,3], col=col4) 

plot(simClo.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim){
  lines(simClo.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,3], col=col5)

dev.off()