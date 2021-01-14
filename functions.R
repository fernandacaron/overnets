
simRangeOver <- function(meanlog, S, truncate=TRUE) {
  library(moments)
  S <- S
  range.sizes <- rlnorm(S, meanlog = meanlog)
  if(truncate==TRUE) range.sizes[range.sizes > 1] <- 1
  rho<-round(mean(range.sizes), digits = 3)
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


simRangeOver2<-function(meanlog, sdlog, S, truncate=TRUE, D.min=D.min, D.max=D.max) {
  S <- S
  range.sizes <- rlnorm(S, meanlog = meanlog, sdlog = sdlog)
  if(truncate==TRUE) range.sizes[range.sizes > D.max] <- D.max
  rho<-mean(log(range.sizes))/(D.max-D.min)
  range.limits <- matrix(ncol = 2, nrow = S)
  midpoint <- runif(S, min = D.min, max = D.max)
  for (i in 1:S) {
    half.rangesizes <- range.sizes/2
    range.limits[i, 1] <- midpoint[i] - half.rangesizes[i]
    range.limits[i, 2] <- midpoint[i] + half.rangesizes[i]
  }
  if(truncate==TRUE) range.limits[range.limits > D.max] <- D.max
  if(truncate==TRUE) range.limits[range.limits < D.min] <- D.min
  
  range.limit.max<-apply(range.limits, 1, max)
  range.limit.min<-apply(range.limits, 1, min)
  overlap<-matrix(nrow = S, ncol = S)
  for(i in 1:S){
    pair.max<-mapply(min, range.limit.max[i], range.limit.max)
    pair.min<-mapply(max, range.limit.min[i], range.limit.min)
    overlap[i,]<-pair.max-pair.min
  }
  overlap[overlap<0]<-0

  res<-list()
  res$rho<-rho
  res$over<-overlap
  return(res)
}

overlap.matrix <- function(range.limits){
  range.limit.max<-apply(range.limits, 1, max)
  range.limit.min<-apply(range.limits, 1, min)
  overlap<-matrix(nrow = nrow(range.limits), ncol = nrow(range.limits))
  for(i in 1:nrow(range.limits)){
    pair.max<-mapply(min, range.limit.max[i], range.limit.max)
    pair.min<-mapply(max, range.limit.min[i], range.limit.min)
    overlap[i,]<-pair.max-pair.min
  }
  overlap[overlap<0]<-0
  return(overlap)
}

overcalc <- function(v1, v2) {
  ov <- min(max(v1), max(v2)) - max(min(v1), min(v2))
  ifelse(ov > 0, ov, 0)
}

calcMetrics <- function(matrix){
  library(sna)
  
  deg<-degree(matrix, gmode = "graph")
  bet<-betweenness(matrix, gmode = "graph")
  clo<-closeness(matrix, gmode = "graph")
  
  res<-list()
  res$degree<-deg
  res$betweeness<-bet
  res$closeness<-clo
  
  return(res)
}
