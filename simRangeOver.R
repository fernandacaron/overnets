# This function generates a matrix of range overlap among a set of S species
# placed along a single dimension with length D. The range of each 
# species was simulated using a lognormal distribution, with a mean of 0 and a sd 
# on a log scale of 1, which generates an average range size of ~1.6. If the 
# simulated range size was larger than D, it was truncated to the value of D. 

simRangeOver <- function(meanlog, S, truncate=TRUE) {
  library(moments)
  S <- S
  range.sizes <- rlnorm(S, meanlog = meanlog)
  rho<-round(mean(range.sizes), digits = 3)
  if(truncate==TRUE) range.sizes[range.sizes > 1] <- 1
  range.limits <- matrix(ncol = 2, nrow = S)
  midpoint <- runif(S, min = 0, max = 1)
  for (i in 1:S) {
    half.rangesizes <- range.sizes/2
    range.limits[i, 1] <- midpoint[i] - half.rangesizes[i]
    range.limits[i, 2] <- midpoint[i] + half.rangesizes[i]
  }
  if(truncate==TRUE) range.limits[range.limits > 1] <- 1
  if(truncate==TRUE) range.limits[range.limits < 0] <- 0
  
  overlap <- matrix(nrow = S, ncol = S)
  overcalc <- function(v1, v2) {
    ov <- min(max(v1), max(v2)) - max(min(v1), min(v2))
    ifelse(ov > 0, ov, 0)
  }
  
  for (i in 1:S) {
    for (j in 1:S) {
      pair <- rbind(range.limits[i, ], range.limits[j, ])
      overlap[i, j] <- overcalc(pair[1, ], pair[2, ])
    }
  }
  res<-list()
  res$rho<-rho
  res$over<-overlap
  return(res)
}
