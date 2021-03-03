#Simulating range overlap distribution
#
#This function generates a matrix of range overlap among a set of S species
#placed along a single dimension with length D (Dmax-Dmin). The range of each 
#species was simulated using a lognormal distribution. If the simulated range 
#size was larger than D, it was truncated to the value of D (if truncate=TRUE). 
#
#meanlog,sdlog: Mean and standard deviation of the distribution on the log scale 
#for the range size simulation.
#S: Number of species to simulate.
#truncation: logical; if TRUE, range sizes will be truncated when outside the Domain.
#Dmin,Dmax: Minimum and maximum size of the Domain.
#
#Returns: A matrix of the range overlap between all species.
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

#Calculating a range overlap matrix
#
#This function generates a matrix of range overlap from a given set of range
#limits.
#
#range.limits: Data frame with the lower and upper limit of the range size in 
#each column, with rows representing the species.
#
#Return: A matrix of the range overlap between all species.
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

#Calculating network metrics of a range overlap matrix
#
#This function calculates network metrics (degree, betweenness, and closeness) 
#of a simulated range overlap matrix for a set of values of mean or standard 
#deviation of the distribution on the log scale for the range size simulation,
#The resulting metrics for a given mean or standard deviation are the average 
#value over 100 repetitions. This function does not work if given mean and 
#standard deviation simultaneously.
#
#meanlog,sdlog: Value or set of values of mean or standard deviation of the 
#distribution on the log scale for the range size simulation
#S: Number of species to be simulated.
#truncation: logical; if TRUE, range sizes will be truncated when outside the Domain
#Return: A list with the degree, betweenness, and closeness values of each species, 
#with each line representing simulation for a given value of meanlog or sdlog.
netStats <- function(meanlog, sdlog, S, Dmin, Dmax, truncate = TRUE) {
	nStats <- function(meanlog, sdlog, S, Dmin, Dmax, truncate = TRUE) {
		nreps <- 100
		res <- list()
		bt <- cl <- dg  <- matrix(ncol = S, nrow = nreps)
		for (i in 1:nreps) {
			x <- simRangeOver(meanlog=meanlog, sdlog=sdlog, S=S, Dmin=Dmin, 
								Dmax=Dmax, truncate = truncate)$overlap
			diag(x)<-NA
			x[x > 0] <- 1
			bt[i, ] <- sort(betweenness(x, gmode = "graph"), decreasing = TRUE)
			cl[i, ] <- sort(closeness(x, gmode = "graph"), decreasing = TRUE)
			dg[i, ] <- sort(degree(x, gmode = "graph"), decreasing = TRUE)
		}
		res$bt <- colMeans(bt)
		res$cl <- colMeans(cl)
		res$dg <- colMeans(dg)
		res
	}
	nreps <- 100

	if(missing(sdlog)){
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:length(meanlog)) {
		x <- nStats(meanlog[i], sdlog=1, S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	} else { if(missing(meanlog)) {
	res <- list()
	res$bt <- res$cl <- res$dg <- matrix(ncol = S, nrow = nreps)
	for (i in 1:length(sdlog)) {
		x <- nStats(meanlog=0, sdlog=sdlog[i], S, Dmin, Dmax, truncate = TRUE)
		res$bt[i, ] <- x$bt
		res$cl[i, ] <- x$cl
		res$dg[i, ] <- x$dg
	}
	res
	}
	}
}