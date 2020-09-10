
genOver <- function(meanlog) {
# This function generates an illustration of the range sizes of different species along
# a single dimension D between 0 and 1, as we vary meanlog
  S <- 50
  range.sizes <- rlnorm(S, meanlog=meanlog)
  rho<-round(mean(range.sizes), digits = 2)
  range.sizes[range.sizes>1]<-1
  range.limits <- matrix(ncol = 2, nrow = S)
  midpoint <- runif(S, min = 0, max = 1)
  for (i in 1:S) {
    half.rangesizes <- range.sizes / 2
    range.limits[i, 1] <- midpoint[i] - half.rangesizes[i]
    range.limits[i, 2] <- midpoint[i] + half.rangesizes[i]
  }
    range.limits[range.limits<0]<-0
    range.limits[range.limits>1]<-1

  plot(
    1,
    xlim = c(0, 1),
    ylim = c(1, S),
    type = "n",
    xlab = "D",
    ylab = "Species ID",
    main = ""

  )
  for (j in 1:S) {
    lines(range.limits[j,], c(j, j))
  }
}
