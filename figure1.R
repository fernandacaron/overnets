rm(list=ls())

#Creating range sizes for 20 species at random, plotting the range sizes 
#in order in the middle of the Domain and then spread across the Domain 
#to show how truncation works.

S<-20
x<-sample(c(0.05, 0.1, 0.2, 0.2, 0.25, 0.3, 0.4, 0.35, 0.6, 0.7, 0.05, 0.1, 0.2, 0.2, 0.25, 0.3, 0.4, 0.35, 0.6, 0.7))
y<-1:20


pdf("Figure1.pdf", width=3)
layout(matrix(1:3, ncol=1))

#FigA
plot(x, ylim=c(0,20), xlim=c(-0.1, 1.1), type="n", xlab="D", ylab="Species")
for(i in 1:20) lines(c(y[i],y[i])~c(0.5-x[i]/2, 0.5+x[i]/2))
text(1,20, "      A", font=3)

#FigB
plot(x, ylim=c(0,20), xlim=c(-0.1, 1.1), type="n", xlab="D", ylab="Species")
for(i in 1:20) {
	xx<-c(0.5-x[i]/2, 0.5+x[i]/2)+runif(1, -0.5, 0.5)
	lines(c(y[i],y[i])~xx)
	abline(v=0, col="red")
	abline(v=1, col="red")
}
text(1,20, "      B", font=3)
#FigC
plot(x, ylim=c(0,20), xlim=c(-0.1, 1.1), type="n", xlab="D", ylab="Species")
for(i in 1:20) {
	xx<-c(0.5-x[i]/2, 0.5+x[i]/2)+runif(1, -0.5, 0.5)
	xx[xx<0]<-0
	xx[xx>1]<-1
	lines(c(y[i],y[i])~xx)
	abline(v=0, col="red")
	abline(v=1, col="red")
}
text(1,20, "      C", font=3)
dev.off()