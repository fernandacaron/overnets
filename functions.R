#' Simulating range overlap distribution
#'
#' This function generates a matrix of range overlap among a set of S species
#' placed along a single dimension with length D (Dmax-Dmin). The range of each 
#' species was simulated using a lognormal distribution. If the simulated range 
#' size was larger than D, it was truncated to the value of D (if truncate=TRUE). 
#'
#' @param meanlog, sdlog Mean and standard deviation of the distribution on the log scale for the range size simulation
#' @param S Number of species to simulate
#' @param truncation logical; if TRUE, range sizes will be truncated when outside the Domain
#' @param Dmin, Dmax Minimum and maximum size of the Domain
#' @return A matrix of the range overlap between all species
#' @export
simRangeOver <- function(meanlog, sdlog, S, truncate=TRUE, Dmin = Dmin, Dmax = Dmax) {
  
  S <- S
  range.sizes <- rlnorm(S, meanlog = meanlog, sdlog = sdlog)
  rho <- (mean(log(range.sizes)))/(log(Dmax - Dmin))
  range.limits <- matrix(ncol = 2, nrow = S)
  midpoint <- runif(S, min = Dmin, max = Dmax)
  for (i in 1:S) {
    half.rangesizes <- range.sizes/2
    range.limits[i, 1] <- midpoint[i] - half.rangesizes[i]
    range.limits[i, 2] <- midpoint[i] + half.rangesizes[i]
  }
  if(truncate==TRUE) range.limits[range.limits > Dmax] <- Dmax
  if(truncate==TRUE) range.limits[range.limits < Dmin] <- Dmin
  
  range.limit.max <- apply(range.limits, 1, max)
  range.limit.min <- apply(range.limits, 1, min)
  overlap <- matrix(nrow = S, ncol = S)
  for(i in 1:S){
    pair.max <- mapply(min, range.limit.max[i], range.limit.max)
    pair.min <- mapply(max, range.limit.min[i], range.limit.min)
    overlap[i,] <- pair.max-pair.min
  }
  overlap[overlap < 0] <- 0
  
  res <- list()
  res$range.sizes <- range.sizes
  res$rho <- rho
  res$overlap <- overlap
  return(res)
}

#' Calculating a range overlap matrix
#'
#' This function generates a matrix of range overlap from a given set of range
#' limits.
#'
#' @param range.limits Lower and higher limit of the range distribution
#' @return A matrix of the range overlap between all species
#' @export
overlap.matrix <- function(range.limits){
  
  range.limit.max <- apply(range.limits, 1, max)
  range.limit.min <- apply(range.limits, 1, min)
  overlap <- matrix(nrow = nrow(range.limits), ncol = nrow(range.limits))
  for(i in 1:nrow(range.limits)){
    pair.max <- mapply(min, range.limit.max[i], range.limit.max)
    pair.min <- mapply(max, range.limit.min[i], range.limit.min)
    overlap[i,] <- pair.max-pair.min
  }
  overlap[overlap < 0] <- 0
  
  return(overlap)
}

  #This function calculates network metrics (degree, betweenness and closeness)
  #of a range overlap matrix. 

#' Calculating network metrics of a range overlap matrix
#'
#' This function calculates network metrics (degree, betweenness and closeness)
#' of a simulated range overlap matrix for a set of values of mean or standard 
#' deviation of the distribution on the log scale for the range size simulation.
#' The resulting metrics for a given mean or standard deviation are the average 
#' value over 100 repetitions.
#'
#' @param fixed.sdlog logical; if TRUE, the range overlap matrices will be simulated 
#' over a set of mean of the distribution on the log scale, with fixed standard deviation; if FALSE,
#' simulation will be performed over a set of standard deviation values, with a fixed mean.
#' @param S Number of species to be simulated
#' @param truncation logical; if TRUE, range sizes will be truncated when outside the Domain
#' @return A list with degree, closenness and betweenness values of each species
#' @export
netStats <- function(fixed.sdlog = TRUE, S, truncate = TRUE) {
	nStats <- function(meanlog, sdlog, S, Dmin, Dmax, truncate = TRUE) {
		nreps <- 100
		res <- list()
		bt <- cl <- dg  <- matrix(ncol = S, nrow = nreps)
		for (i in 1:nreps) {
			x <- simRangeOver(meanlog=meanlog, sdlog=sdlog, S=S, Dmin=Dmin, 
								Dmax=Dmax, truncate = truncate)$overlap
			diag(x)<-NA
			x[x > 0] <- 1
			bt[i, ] <- sort(sna::betweenness(x, gmode = "graph"), decreasing = TRUE)
			cl[i, ] <- sort(sna::closeness(x, gmode = "graph"), decreasing = TRUE)
			dg[i, ] <- sort(sna::degree(x, gmode = "graph"), decreasing = TRUE)
		}
		res$bt <- colMeans(bt)
		res$cl <- colMeans(cl)
		res$dg <- colMeans(dg)
		res
	}
	nreps <- 100
	meanlog <- seq(from = -4, to = 0.5, by = 0.045)
	sdlog<-seq(from = 0, to = 10, by = 0.1)
	Dmin <- 0
	Dmax <- exp(1)
	if(fixed.sdlog==TRUE){
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:nreps) {
	print(i)
		x <- nStats(meanlog[i], sdlog=1, S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	} else {
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:nreps) {
	print(i)
		x <- nStats(meanlog=0, sdlog=sdlog[i], S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	}
}