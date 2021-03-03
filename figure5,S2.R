rm(list=ls())

library(MASS)
library(sna)

source("functions.R")

#This code will calculate range size for an empirical dataset, estimate the
#parameters of a fitted lognormal distribution to the range sizes, calculate
#the network metrics and simulate range overlap matrices to calculate the 
#expected values

dat<-read.csv("adata.csv")

#This function will subset the data and calculate the range sizes for a given 
#location 
subsetData<-function(data, location){
  dat<-cbind(data[colnames(data)==location], data[,c(3,4,10,12)])
  dat<-dat[dat[,1]=="Y",]
  rownames(dat)<-paste(dat[,2], dat[,3], sep="_")
  dat<-dat[,c(4,5)]
  dat[,c(1:2)]<-apply(dat[,c(1:2)], 2, function(x) as.numeric((x)))
  dat<-na.omit(dat)
  
  range.sizes<-numeric()
  for(i in 1:nrow(dat)){
    range.sizes[i]<-dat$MAX[i]-dat$MIN[i]
  }
  names(range.sizes)<-rownames(dat)
  range.sizes<-range.sizes[range.sizes>0]
  
  dat<-dat[rownames(dat) %in% names(range.sizes),]  
  
  subset<-list()
  subset$range.sizes<-range.sizes
  subset$dat<-dat
  return(subset)
}

datMAH<-subsetData(dat, "MAH")$dat
datCDH<-subsetData(dat, "CDH")$dat
datCAN<-subsetData(dat, "CAN")$dat
datNAN<-subsetData(dat, "NAN")$dat
datTEP<-subsetData(dat, "TEP")$dat

range.sizes.MAH<-subsetData(dat, "MAH")$range.sizes
range.sizes.CDH<-subsetData(dat, "CDH")$range.sizes
range.sizes.CAN<-subsetData(dat, "CAN")$range.sizes
range.sizes.NAN<-subsetData(dat, "NAN")$range.sizes
range.sizes.TEP<-subsetData(dat, "TEP")$range.sizes

#Taking the number of species of each location and the rho values
S.MAH<-length(range.sizes.MAH)
S.CDH<-length(range.sizes.CDH)
S.CAN<-length(range.sizes.CAN)
S.NAN<-length(range.sizes.NAN)
S.TEP<-length(range.sizes.TEP)

rho.MAH<-mean(log(range.sizes.MAH))/log((max(datMAH$MAX)-min(datMAH$MIN)))
rho.CDH<-mean(log(range.sizes.CDH))/log((max(datCDH$MAX)-min(datCDH$MIN)))
rho.CAN<-mean(log(range.sizes.CAN))/log((max(datCAN$MAX)-min(datCAN$MIN)))
rho.NAN<-mean(log(range.sizes.NAN))/log((max(datNAN$MAX)-min(datNAN$MIN)))
rho.TEP<-mean(log(range.sizes.TEP))/log((max(datTEP$MAX)-min(datTEP$MIN)))

#Fitting a lognormal distribution to the data and using the estimates to 
#simulate the range overlap matrices
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

#Calculating the overlap matrix
overlap.MAH<-overlap.matrix(datMAH)
overlap.CDH<-overlap.matrix(datCDH)
overlap.CAN<-overlap.matrix(datCAN)
overlap.NAN<-overlap.matrix(datNAN)
overlap.TEP<-overlap.matrix(datTEP)

#Making the overlap matrix binary for calculating the network metrics
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

#Calculating the average number of overlapping species in the datasets in the 
#same binary matrices
meanNumOverSpp.MAH<-mean(rowSums(over.MAH))
meanNumOverSpp.CDH<-mean(rowSums(over.CDH))
meanNumOverSpp.CAN<-mean(rowSums(over.CAN))
meanNumOverSpp.NAN<-mean(rowSums(over.NAN))
meanNumOverSpp.TEP<-mean(rowSums(over.TEP))

#Getting only the range overlaps higher than 0 of the complete matrices to see
#the shape of the distribution
overlap.MAH<-overlap.MAH[overlap.MAH>0]
overlap.CDH<-overlap.CDH[overlap.CDH>0]
overlap.CAN<-overlap.CAN[overlap.CAN>0]
overlap.NAN<-overlap.NAN[overlap.NAN>0]
overlap.TEP<-overlap.TEP[overlap.TEP>0]

#Simulations based on the empirical datasets

#Setting the number o simulations; nsim1 corresponds to the number of range 
#overlap matrices and nsim2 corresponds to how many simulated matrices will
#get your network metrics calculated
nsim1<-100
nsim2<-50

#This function will take all the parameters of each location and perform the 
#simulations of the average number of overlapping species, rho values, the 
#variance of range overlap and variance of range size, the mean range overlap,
#and the metricsof the networks.
sim<-function(meanlog, sdlog, S, truncate=TRUE, Dmin, Dmax, nsim1, nsim2){
  simMeanNumOverSpp<-simRho<-simVarOver<-numeric()
  simOver<-simMeanOver<-list()
  for(i in 1:nsim1){
    dat<-simRangeOver(meanlog=meanlog, sdlog=sdlog, S=S, truncate=truncate, 
                      Dmin=Dmin, Dmax=Dmax)
    overlap<-dat$overlap
    over<-dat$overlap
    diag(overlap)<-0
    overlap[overlap>0]<-1
    simMeanNumOverSpp[i]<-mean(rowSums(overlap))
    simRho[i]<-dat$rho
    
    simOver[[i]]<-overlap
    
    diag(over)<-0
    over<-over[over>0]
    simMeanOver[[i]]<-log(over)
    simVarOver[i]<-var(as.vector(log(over)), na.rm = TRUE)
  }
  simMeanOver<-do.call(c, simMeanOver)
  
  simDeg<-simBet<-simClo<-matrix(nrow=S, ncol=nsim2)
  for(i in 1:nsim2){
    simDeg[,i]<-sort(degree(simOver[[i]], gmode="graph"), decreasing=T)
    simBet[,i]<-sort(betweenness(simOver[[i]], gmode="graph"), decreasing=T)
    simClo[,i]<-sort(closeness(simOver[[i]], gmode="graph"), decreasing=T)
  }
  
  sim<-list()
  
  sim$MeanNumOverSpp<-simMeanNumOverSpp
  sim$Rho<-simRho
  sim$VarOver<-simVarOver
  sim$MeanOver<-simMeanOver
  sim$Deg<-simDeg
  sim$Bet<-simBet
  sim$Clo<-simClo
  
  return(sim)
}

sim.MAH.T<-sim(meanlog=meanlog.MAH, sdlog=sdlog.MAH, S=S.MAH, truncate=FALSE, 
               Dmin=min(datMAH$MIN), Dmax=max(datMAH$MAX), nsim1, nsim2)
simMeanNumOverSpp.MAH.T<-sim.MAH.T$MeanNumOverSpp
simRho.MAH.T<-sim.MAH.T$Rho
simVarOver.MAH.T<-sim.MAH.T$VarOver
simMeanOver.MAH.T<-sim.MAH.T$MeanOver
simDeg.MAH.T<-sim.MAH.T$Deg
simBet.MAH.T<-sim.MAH.T$Bet
simClo.MAH.T<-sim.MAH.T$Clo

sim.CDH.T<-sim(meanlog=meanlog.CDH, sdlog=sdlog.CDH, S=S.CDH, truncate=FALSE, 
               Dmin=min(datCDH$MIN), Dmax=max(datCDH$MAX), nsim1, nsim2)
simMeanNumOverSpp.CDH.T<-sim.CDH.T$MeanNumOverSpp
simRho.CDH.T<-sim.CDH.T$Rho
simVarOver.CDH.T<-sim.CDH.T$VarOver
simMeanOver.CDH.T<-sim.CDH.T$MeanOver
simDeg.CDH.T<-sim.CDH.T$Deg
simBet.CDH.T<-sim.CDH.T$Bet
simClo.CDH.T<-sim.CDH.T$Clo

sim.CAN.T<-sim(meanlog=meanlog.CAN, sdlog=sdlog.CAN, S=S.CAN, truncate=FALSE, 
               Dmin=min(datCAN$MIN), Dmax=max(datCAN$MAX), nsim1, nsim2)
simMeanNumOverSpp.CAN.T<-sim.CAN.T$MeanNumOverSpp
simRho.CAN.T<-sim.CAN.T$Rho
simVarOver.CAN.T<-sim.CAN.T$VarOver
simMeanOver.CAN.T<-sim.CAN.T$MeanOver
simDeg.CAN.T<-sim.CAN.T$Deg
simBet.CAN.T<-sim.CAN.T$Bet
simClo.CAN.T<-sim.CAN.T$Clo

sim.NAN.T<-sim(meanlog=meanlog.NAN, sdlog=sdlog.NAN, S=S.NAN, truncate=FALSE, 
               Dmin=min(datNAN$MIN), Dmax=max(datNAN$MAX), nsim1, nsim2)
simMeanNumOverSpp.NAN.T<-sim.NAN.T$MeanNumOverSpp
simRho.NAN.T<-sim.NAN.T$Rho
simVarOver.NAN.T<-sim.NAN.T$VarOver
simMeanOver.NAN.T<-sim.NAN.T$MeanOver
simDeg.NAN.T<-sim.NAN.T$Deg
simBet.NAN.T<-sim.NAN.T$Bet
simClo.NAN.T<-sim.NAN.T$Clo

sim.TEP.T<-sim(meanlog=meanlog.TEP, sdlog=sdlog.TEP, S=S.TEP, truncate=FALSE, 
               Dmin=min(datTEP$MIN), Dmax=max(datTEP$MAX), nsim1, nsim2)
simMeanNumOverSpp.TEP.T<-sim.TEP.T$MeanNumOverSpp
simRho.TEP.T<-sim.TEP.T$Rho
simVarOver.TEP.T<-sim.TEP.T$VarOver
simMeanOver.TEP.T<-sim.TEP.T$MeanOver
simDeg.TEP.T<-sim.TEP.T$Deg
simBet.TEP.T<-sim.TEP.T$Bet
simClo.TEP.T<-sim.TEP.T$Clo

pdf("Figure5.pdf", width=12, height=14)
layout(matrix(1:35, nrow=7, ncol=5, byrow=T))
par(mar=c(4, 4, 1, 1))

col1<-rgb(33/255,12/255,74/255, 0.8)
col2<-rgb(137/255,34/255,106/255, 0.8)
col3<-rgb(187/255,55/255,84/255, 0.8)
col4<-rgb(227/255,89/255,50/255, 0.8)
col5<-rgb(249/255,140/255,10/255, 0.8)

#First row: Range sizes distributed in the Domain
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

#Second row: Range size distribution
hist(log(range.sizes.MAH), col=col1, border=col1, 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.CDH), col=col2, border=col2, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.CAN), col=col3, border=col3, ylab="",  
     xlim=c(0,10), xlab="Range size", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.NAN), col=col4, border=col4, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.TEP), col=col5, border=col5, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

#Third row: Observed and expected range overlap distribution with the respective
#means
hist(log(overlap.MAH), col=col1, border=col1, main="", xlab="", ylab="Density",
     breaks=seq(from=-10, to=10, by=0.2), freq=F, xlim=c(0,10), ylim=c(0,1))
abline(v=mean(log(overlap.MAH), na.rm=T), col=col1)
hist(simMeanOver.MAH.NT, col=rgb(0.7,0.7,0.7,0.4), freq=F,
     border=rgb(0.7,0.7,0.7,0.4), add=T, 
     breaks=seq(from=min(simMeanOver.MAH.NT)-2, to=max(simMeanOver.MAH.NT)+1, 
                by=0.2))
abline(v=mean(simMeanOver.MAH.NT), col=rgb(0.7,0.7,0.7))

hist(log(overlap.CDH), col=col2, border=col2, main="", xlab="", ylab="", 
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.CDH), na.rm=T), col=col2)
hist(simMeanOver.CDH.NT, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), 
     breaks=seq(from=min(simMeanOver.CDH.NT)-2, to=max(simMeanOver.CDH.NT)+2, 
                by=0.2))
abline(v=mean(simMeanOver.CDH.NT), col=rgb(0.7,0.7,0.7))

hist(log(overlap.CAN), col=col3, border=col3, main="", xlab="Range overlap", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.CAN), na.rm=T), col=col3)
hist(simMeanOver.CAN.NT, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), 
     breaks=seq(from=min(simMeanOver.CAN.NT)-2, to=max(simMeanOver.CAN.NT)+2, 
                by=0.2))
abline(v=mean(simMeanOver.CAN.NT), col=rgb(0.7,0.7,0.7))

hist(log(overlap.NAN), col=col4, border=col4, main="", xlab="", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.NAN), na.rm=T), col=col4)
hist(simMeanOver.NAN.NT, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), breaks=seq(from=min(simMeanOver.NAN.NT)-2, 
                                             to=max(simMeanOver.NAN.NT)+2, by=0.2))
abline(v=mean(simMeanOver.NAN.NT), col=rgb(0.7,0.7,0.7))

hist(log(overlap.TEP), col=col5, border=col5, main="", xlab="", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.TEP), na.rm=T), col=col5)
hist(simMeanOver.TEP.NT, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F, 
     border=rgb(0.7,0.7,0.7,0.4), breaks=seq(from=min(simMeanOver.TEP.NT)-2, 
                                             to=max(simMeanOver.TEP.NT)+2, by=0.2))
abline(v=mean(simMeanOver.TEP.NT), col=rgb(0.7,0.7,0.7))

#Fourth row: Expected distribution of the average number of overlapping species 
#and the observed line
hist(simMeanNumOverSpp.MAH.NT, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), xlab="", main="", xlim=c(200,400),
     breaks=seq(from=200, to=400, by=3.5))
abline(v=meanNumOverSpp.MAH, col=col1)

hist(simMeanNumOverSpp.CDH.NT, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(140,270),
     breaks=seq(from=140, to=270, by=2))
abline(v=meanNumOverSpp.CDH, col=col2)

hist(simMeanNumOverSpp.CAN.NT, col=rgb(0.7,0.7,0.7,0.6), xlim=c(350,600),
     border=rgb(0.7,0.7,0.7,0.6), ylab="", main="",
     xlab="Average number of overlapping species", 
     breaks=seq(from=350, to=600, by=4))
abline(v=meanNumOverSpp.CAN, col=col3)

hist(simMeanNumOverSpp.NAN.NT, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(350,600),
     breaks=seq(from=350, to=600, by=4))
abline(v=meanNumOverSpp.NAN, col=col4)

hist(simMeanNumOverSpp.TEP.NT, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(60,140),
     breaks=seq(from=60, to=140, by=1.5))
abline(v=meanNumOverSpp.TEP, col=col5)

#Fifth row: Degree distribution
plot(simDeg.MAH.NT[,1], main="", ylab="Degree", xlab="", type="l", 
     ylim=range(simDeg.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,1], col=col1) 

plot(simDeg.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,1], col=col2) 

plot(simDeg.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simDeg.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,1], col=col3) 

plot(simDeg.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,1], col=col4) 

plot(simDeg.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,1], col=col5)

#Sixth row: Betweenness distribution.
plot(simBet.MAH.NT[,1], main="", ylab="Betweenness", xlab="", type="l", 
     ylim=range(simBet.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,2], col=col1) 

plot(simBet.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,2], col=col2) 

plot(simBet.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simBet.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,2], col=col3) 

plot(simBet.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,2], col=col4) 

plot(simBet.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,2], col=col5)

#Seventh row: Closeness distribution.
plot(simClo.MAH.NT[,1], main="", ylab="Closeness", xlab="", type="l", 
     ylim=range(simClo.MAH.NT), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.MAH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,3], col=col1) 

plot(simClo.CDH.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.CDH.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.CDH.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,3], col=col2) 

plot(simClo.CAN.NT[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simClo.CAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.CAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,3], col=col3) 

plot(simClo.NAN.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.NAN.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.NAN.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,3], col=col4) 

plot(simClo.TEP.NT[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.TEP.NT[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.TEP.NT[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,3], col=col5)

dev.off()

#Repeting the simulations with truncation

sim.MAH.T<-sim(meanlog=meanlog.MAH, sdlog=sdlog.MAH, S=S.MAH, truncate=TRUE, 
               Dmin=min(datMAH$MIN), Dmax=max(datMAH$MAX), nsim1, nsim2)
simMeanNumOverSpp.MAH.T<-sim.MAH.T$MeanNumOverSpp
simRho.MAH.T<-sim.MAH.T$Rho
simVarOver.MAH.T<-sim.MAH.T$VarOver
simMeanOver.MAH.T<-sim.MAH.T$MeanOver
simDeg.MAH.T<-sim.MAH.T$Deg
simBet.MAH.T<-sim.MAH.T$Bet
simClo.MAH.T<-sim.MAH.T$Clo

sim.CDH.T<-sim(meanlog=meanlog.CDH, sdlog=sdlog.CDH, S=S.CDH, truncate=TRUE, 
               Dmin=min(datCDH$MIN), Dmax=max(datCDH$MAX), nsim1, nsim2)
simMeanNumOverSpp.CDH.T<-sim.CDH.T$MeanNumOverSpp
simRho.CDH.T<-sim.CDH.T$Rho
simVarOver.CDH.T<-sim.CDH.T$VarOver
simMeanOver.CDH.T<-sim.CDH.T$MeanOver
simDeg.CDH.T<-sim.CDH.T$Deg
simBet.CDH.T<-sim.CDH.T$Bet
simClo.CDH.T<-sim.CDH.T$Clo

sim.CAN.T<-sim(meanlog=meanlog.CAN, sdlog=sdlog.CAN, S=S.CAN, truncate=TRUE, 
               Dmin=min(datCAN$MIN), Dmax=max(datCAN$MAX), nsim1, nsim2)
simMeanNumOverSpp.CAN.T<-sim.CAN.T$MeanNumOverSpp
simRho.CAN.T<-sim.CAN.T$Rho
simVarOver.CAN.T<-sim.CAN.T$VarOver
simMeanOver.CAN.T<-sim.CAN.T$MeanOver
simDeg.CAN.T<-sim.CAN.T$Deg
simBet.CAN.T<-sim.CAN.T$Bet
simClo.CAN.T<-sim.CAN.T$Clo

sim.NAN.T<-sim(meanlog=meanlog.NAN, sdlog=sdlog.NAN, S=S.NAN, truncate=TRUE, 
               Dmin=min(datNAN$MIN), Dmax=max(datNAN$MAX), nsim1, nsim2)
simMeanNumOverSpp.NAN.T<-sim.NAN.T$MeanNumOverSpp
simRho.NAN.T<-sim.NAN.T$Rho
simVarOver.NAN.T<-sim.NAN.T$VarOver
simMeanOver.NAN.T<-sim.NAN.T$MeanOver
simDeg.NAN.T<-sim.NAN.T$Deg
simBet.NAN.T<-sim.NAN.T$Bet
simClo.NAN.T<-sim.NAN.T$Clo

sim.TEP.T<-sim(meanlog=meanlog.TEP, sdlog=sdlog.TEP, S=S.TEP, truncate=TRUE, 
               Dmin=min(datTEP$MIN), Dmax=max(datTEP$MAX), nsim1, nsim2)
simMeanNumOverSpp.TEP.T<-sim.TEP.T$MeanNumOverSpp
simRho.TEP.T<-sim.TEP.T$Rho
simVarOver.TEP.T<-sim.TEP.T$VarOver
simMeanOver.TEP.T<-sim.TEP.T$MeanOver
simDeg.TEP.T<-sim.TEP.T$Deg
simBet.TEP.T<-sim.TEP.T$Bet
simClo.TEP.T<-sim.TEP.T$Clo

pdf("figures/FigureS2.pdf", width=12, height=14)
layout(matrix(1:35, nrow=7, ncol=5, byrow=T))
par(mar=c(4, 4, 1, 1))

col1<-rgb(33/255,12/255,74/255, 0.8)
col2<-rgb(137/255,34/255,106/255, 0.8)
col3<-rgb(187/255,55/255,84/255, 0.8)
col4<-rgb(227/255,89/255,50/255, 0.8)
col5<-rgb(249/255,140/255,10/255, 0.8)

#First row: Range sizes distributed in the Domain
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


#Second row: Range size distribution
hist(log(range.sizes.MAH), col=col1, border=col1, 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.CDH), col=col2, border=col2, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.CAN), col=col3, border=col3, ylab="",  
     xlim=c(0,10), xlab="Range size", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.NAN), col=col4, border=col4, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

hist(log(range.sizes.TEP), col=col5, border=col5, ylab="", 
     xlim=c(0,10), xlab="", main="", breaks=seq(from=-10, to=10, by=0.2))

#Third row: Observed and expected range overlap distribution with the respective
#means
hist(log(overlap.MAH), col=col1, border=col1, main="", xlab="", ylab="Density",
     breaks=seq(from=-10, to=10, by=0.2), freq=F, xlim=c(0,10), ylim=c(0,1))
abline(v=mean(log(overlap.MAH), na.rm=T), col=col1)
hist(simMeanOver.MAH.T, col=rgb(0.7,0.7,0.7,0.4), freq=F,
     border=rgb(0.7,0.7,0.7,0.4), add=T, 
     breaks=seq(from=min(simMeanOver.MAH.T)-2, to=max(simMeanOver.MAH.T)+2, 
                by=0.2))
abline(v=mean(simMeanOver.MAH.T), col=rgb(0.7,0.7,0.7))

hist(log(overlap.CDH), col=col2, border=col2, main="", xlab="", ylab="", 
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.CDH), na.rm=T), col=col2)
hist(simMeanOver.CDH.T, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), 
     breaks=seq(from=min(simMeanOver.CDH.T)-2, to=max(simMeanOver.CDH.T)+2, 
                by=0.2))
abline(v=mean(simMeanOver.CDH.T), col=rgb(0.7,0.7,0.7))

hist(log(overlap.CAN), col=col3, border=col3, main="", xlab="Range overlap", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.CAN), na.rm=T), col=col3)
hist(simMeanOver.CAN.T, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), 
     breaks=seq(from=min(simMeanOver.CAN.T)-2, to=max(simMeanOver.CAN.T)+2, 
                by=0.2))
abline(v=mean(simMeanOver.CAN.T), col=rgb(0.7,0.7,0.7))

hist(log(overlap.NAN), col=col4, border=col4, main="", xlab="", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.NAN), na.rm=T), col=col4)
hist(simMeanOver.NAN.T, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F,
     border=rgb(0.7,0.7,0.7,0.4), breaks=seq(from=min(simMeanOver.NAN.T)-2, 
                                             to=max(simMeanOver.NAN.T)+2, by=0.2))
abline(v=mean(simMeanOver.NAN.T), col=rgb(0.7,0.7,0.7))

hist(log(overlap.TEP), col=col5, border=col5, main="", xlab="", ylab="",
     breaks=seq(from=-10, to=10, by=0.2), ylim=c(0,1), freq=F, xlim=c(0,10))
abline(v=mean(log(overlap.TEP), na.rm=T), col=col5)
hist(simMeanOver.TEP.T, col=rgb(0.7,0.7,0.7,0.4), add=T, freq=F, 
     border=rgb(0.7,0.7,0.7,0.4), breaks=seq(from=min(simMeanOver.TEP.T)-2, 
                                             to=max(simMeanOver.TEP.T)+2, by=0.2))
abline(v=mean(simMeanOver.TEP.T), col=rgb(0.7,0.7,0.7))

#Fourth row: Expected distribution of the average number of overlapping species 
#and the observed line
hist(simMeanNumOverSpp.MAH.T, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), xlab="", main="", xlim=c(200,400),
     breaks=seq(from=200, to=400, by=3.5))
abline(v=meanNumOverSpp.MAH, col=col1)

hist(simMeanNumOverSpp.CDH.T, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(140,270),
     breaks=seq(from=140, to=270, by=2))
abline(v=meanNumOverSpp.CDH, col=col2)

hist(simMeanNumOverSpp.CAN.T, col=rgb(0.7,0.7,0.7,0.6), xlim=c(350,600),
     border=rgb(0.7,0.7,0.7,0.6), ylab="", main="",
     xlab="Average number of overlapping species", 
     breaks=seq(from=350, to=600, by=4))
abline(v=meanNumOverSpp.CAN, col=col3)

hist(simMeanNumOverSpp.NAN.T, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(350,600),
     breaks=seq(from=350, to=600, by=4))
abline(v=meanNumOverSpp.NAN, col=col4)

hist(simMeanNumOverSpp.TEP.T, col=rgb(0.7,0.7,0.7,0.6), 
     border=rgb(0.7,0.7,0.7,0.6), ylab="", xlab="", main="", xlim=c(60,140),
     breaks=seq(from=60, to=140, by=1.5))
abline(v=meanNumOverSpp.TEP, col=col5)


#Fifth row: Degree distribution
plot(simDeg.MAH.T[,1], main="", ylab="Degree", xlab="", type="l", 
     ylim=range(simDeg.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,1], col=col1) 

plot(simDeg.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,1], col=col2) 

plot(simDeg.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simDeg.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,1], col=col3) 

plot(simDeg.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,1], col=col4) 

plot(simDeg.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simDeg.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simDeg.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,1], col=col5)

#Sixth row: Betweenness distribution.
plot(simBet.MAH.T[,1], main="", ylab="Betweenness", xlab="", type="l", 
     ylim=range(simBet.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,2], col=col1) 

plot(simBet.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,2], col=col2) 

plot(simBet.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simBet.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,2], col=col3) 

plot(simBet.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,2], col=col4) 

plot(simBet.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simBet.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simBet.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,2], col=col5)

#Seventh row: Closeness distribution.
plot(simClo.MAH.T[,1], main="", ylab="Closeness", xlab="", type="l", 
     ylim=range(simClo.MAH.T), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.MAH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netMAH[,3], col=col1) 

plot(simClo.CDH.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.CDH.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.CDH.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCDH[,3], col=col2) 

plot(simClo.CAN.T[,1], ylab="", xlab="rank", type="l", 
     ylim=range(simClo.CAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.CAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netCAN[,3], col=col3) 

plot(simClo.NAN.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.NAN.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.NAN.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netNAN[,3], col=col4) 

plot(simClo.TEP.T[,1], ylab="", xlab="", type="l", 
     ylim=range(simClo.TEP.T[,1]), col=rgb(0.7,0.7,0.7,0.4)) 
for(i in 2:nsim2){
  lines(simClo.TEP.T[,i], col=rgb(0.7,0.7,0.7,0.4))
}
lines(netTEP[,3], col=col5)

dev.off()