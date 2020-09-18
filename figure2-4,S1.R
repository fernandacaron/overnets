rm(list = ls())
source("fun/genOver.R")
source("fun/simRangeOver.R")

library(moments)

final1<- final2 <- as.data.frame(matrix(ncol = 14, nrow = 100))
colnames(final1) <-colnames(final2) <- c("rho",
"mean.range.size",        "var.range.size",       "sk.range.size",
"mean.log.range.size",    "var.log.range.size",    "sk.log.range.size",
"mean.range.overlap",     "var.range.overlap",     "sk.overlap",
"mean.log.range.overlap",     "var.log.range.overlap", "sk.log.overlap", 
"prop.zeros.overlap")

rho <- seq(from = 0.01, to = 1, by = 0.01)
D <- 1/rho
S<-1000
for (i in 1:100){
x<-simRangeOver(D=D[i], S=S, truncate=TRUE)
overlap<-x

#range size stats
mean.range.size<-mean(diag(x))
var.range.size<-var(diag(x))
sk.range.size<-skewness(diag(x))
mean.log.range.size<-mean(log(diag(x)))
var.log.range.size<-var(log(diag(x)))
sk.log.range.size<-skewness(log(diag(x)))

#range overlap stats
diag(x) <- NA
x <- matrix(x[x > 0], ncol = 1)
mean.range.overlap <- mean(x, na.rm = TRUE)
mean.log.range.overlap <- mean(log(x), na.rm = TRUE)
var.range.overlap <- var(x, na.rm = TRUE)
var.log.range.overlap <- var(log(x), na.rm = TRUE)
sk.overlap <- skewness(x, na.rm = TRUE)
sk.log.overlap <- skewness(log(x), na.rm = TRUE)
prop.zeros.overlap <- table(overlap[upper.tri(overlap)])[1]/length(overlap[upper.tri(overlap)])

res<-cbind(rho[i],
mean.range.size,        var.range.size,       sk.range.size,
mean.log.range.size,    var.log.range.size,    sk.log.range.size,
mean.range.overlap,     var.range.overlap,     sk.overlap,
mean.log.range.overlap, var.log.range.overlap, sk.log.overlap, 
prop.zeros.overlap)

final1[i,]<-res
}
for (i in 1:100){
x<-simRangeOver(D=D[i], S=S, truncate=FALSE)
overlap<-x

#range size stats
mean.range.size<-mean(diag(x))
var.range.size<-var(diag(x))
sk.range.size<-skewness(diag(x))
mean.log.range.size<-mean(log(diag(x)))
var.log.range.size<-var(log(diag(x)))
sk.log.range.size<-skewness(log(diag(x)))

#range overlap stats
diag(x) <- NA
x <- matrix(x[x > 0], ncol = 1)
mean.range.overlap <- mean(x, na.rm = TRUE)
mean.log.range.overlap <- mean(log(x), na.rm = TRUE)
var.range.overlap <- var(x, na.rm = TRUE)
var.log.range.overlap <- var(log(x), na.rm = TRUE)
sk.overlap <- skewness(x, na.rm = TRUE)
sk.log.overlap <- skewness(log(x), na.rm = TRUE)
prop.zeros.overlap <- table(overlap[upper.tri(overlap)])[1]/length(overlap[upper.tri(overlap)])

res<-cbind(rho[i],
mean.range.size,        var.range.size,       sk.range.size,
mean.log.range.size,    var.log.range.size,    sk.log.range.size,
mean.range.overlap,     var.range.overlap,     sk.overlap,
mean.log.range.overlap, var.log.range.overlap, sk.log.overlap, 
prop.zeros.overlap)

final2[i,]<-res
}


pdf("figures/Figure4.pdf", width=5, height=5)
plot(log(prop.zeros.overlap)~rho, data=final1, xlab = expression(rho), ylab = "log of proportion of isolated species", col = rgb(0, 
	0, 139, 120, maxColorValue = 255), pch = 16)
points(log(prop.zeros.overlap)~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
legend(0.6, 0.9, pch=c(16,16), col=c(rgb(0, 0, 139, 120, maxColorValue = 255), rgb(239, 106, 80, 120, maxColorValue = 255)), c("Truncated", "Non-truncated"), cex=.8)
dev.off()


pdf("figures/Figure2.pdf", width=5, height=8)
layout(matrix(1:6, ncol=2))
par(mar = c(4, 4, 2, 4))
plot(mean.log.range.size ~rho, data=final1, xlab = expression(rho), main="Range size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.log.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(var.log.range.size ~rho, data=final1, xlab = expression(rho), ylab = "Variance", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(var.log.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(sk.log.range.size ~rho, data=final1, xlab = expression(rho), ylab = "Skewness", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(sk.log.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(mean.log.range.overlap ~rho, data=final1, xlab = expression(rho), main="Range overlap size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.log.range.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(var.log.range.overlap ~rho, data=final1, xlab = expression(rho), ylab = "Variance", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(var.log.range.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(sk.log.overlap ~rho, data=final1, xlab = expression(rho), ylab = "Skewness", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(sk.log.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
dev.off()


pdf("figures/FigureS1.pdf", width=5, height=8)
layout(matrix(1:6, ncol=2))
par(mar = c(4, 4, 2, 4))
plot(mean.range.size ~rho, data=final1, xlab = expression(rho), main="Range size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(var.range.size ~rho, data=final1, xlab = expression(rho), ylab = "Variance", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(var.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(sk.range.size ~rho, data=final1, xlab = expression(rho), ylab = "Skewness", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(sk.range.size ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(mean.range.overlap ~rho, data=final1, xlab = expression(rho), main="Range overlap size", ylab = "Mean", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(mean.range.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(var.range.overlap ~rho, data=final1, xlab = expression(rho), ylab = "Variance", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(var.range.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)

plot(sk.overlap ~rho, data=final1, xlab = expression(rho), ylab = "Skewness", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
points(sk.overlap ~rho, data=final2, col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
dev.off()



regEq <- function(lmObj, dig) {
    gsub(":", "*", 
        paste0(
            names(lmObj$model)[1]," = ",
            paste0(
                c(round(lmObj$coef[1], dig), round(sign(lmObj$coef[-1])*lmObj$coef[-1], dig)),
                c("", rep("*", length(lmObj$coef)-1)),
                paste0(c("", names(lmObj$coef)[-1]), c(ifelse(sign(lmObj$coef)[-1]==1," + "," - "), "")),
                collapse=""
            )
        )
    )
}

pdf("figures/Figure3.pdf", width=5, height=8)
layout(matrix(1:6, ncol=2, byrow=TRUE))
par(mar = c(4, 4, 2, 4))

plot(mean.log.range.overlap ~ mean.log.range.size, data=final1, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.log.range.overlap ~ mean.log.range.size, data=final1))
plot(mean.log.range.overlap ~ mean.log.range.size, data=final2, xlab="Range size", main="", ylab = "Range overlap size", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
abline(lm(mean.log.range.overlap ~ mean.log.range.size, data=final2))


plot(var.log.range.overlap~var.log.range.size,data=final1,xlab="Variance in range size",main="",ylab="Variance in range overlap", col=rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(var.log.range.overlap~var.log.range.size, data=final1))
plot(var.log.range.overlap~var.log.range.size,data=final2,xlab="Variance in range size",main="",ylab="Variance in range overlap", col=rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
#abline(lm(var.log.range.overlap~var.log.range.size, data=final2))


plot(sk.log.overlap ~ sk.log.range.size, data=final1, xlab = "Skewness in range overlap", main="", ylab = "Skewness in range overlap", col = rgb(0, 0, 139, 120, maxColorValue = 255), pch = 16)
abline(lm(sk.log.overlap ~ sk.log.range.size, data=final1))
plot(sk.log.overlap ~ sk.log.range.size, data=final2, xlab = "Skewness in range overlap", main="", ylab = "Skewness in range overlap", col = rgb(239, 106, 80, 120, maxColorValue = 255), pch = 16)
#abline(lm(sk.log.overlap ~ sk.log.range.size, data=final2))
dev.off()


m1<-lm(mean.log.range.overlap ~ mean.log.range.size, data=final1)
m2<-lm(mean.log.range.overlap ~ mean.log.range.size, data=final2)
m3<-lm(var.log.range.overlap~var.log.range.size,data=final1)
m4<-lm(var.log.range.overlap~var.log.range.size,data=final2)
m5<-lm(sk.log.overlap ~ sk.log.range.size, data=final1)
m6<-lm(sk.log.overlap ~ sk.log.range.size, data=final2)

eqs<-rbind(
c(regEq(m1, 2), round(summary(m1)$r.squared, 2)), 
c(regEq(m2, 2), round(summary(m2)$r.squared, 2)),
c(regEq(m3, 2), round(summary(m3)$r.squared, 2)),
c(regEq(m4, 2), round(summary(m4)$r.squared, 2)),
c(regEq(m5, 2), round(summary(m5)$r.squared, 2)),
c(regEq(m6, 2), round(summary(m6)$r.squared, 2))
)
