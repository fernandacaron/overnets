rm(list = ls())
library(rgeos)
library(maptools)
library(sna)

overCalc <- function() {
	#This function will take all of the shapefiles on the working directory and calculates both the
	#incidence matrix and the quantitative range overlap matrix. Please notice that only shapefiles
	#should be on the directory
	
	spp <- grep("\\.shp$", list.files(recursive = T), value = T)
	spp <- substr(spp, 1, nchar(spp) - 4)

	#Creating incidence matrix
	overBin <- matrix(ncol = length(spp), nrow = length(spp))
	rownames(overBin) <- colnames(overBin) <- spp

	for (i in 1:length(spp)) {
		for (j in 1:length(spp)) {
			sh_i <- readShapePoly(spp[i])
			sh_j <- readShapePoly(spp[j])
			ifelse(gIntersects(sh_i, sh_j), overBin[i, j] <- 1, overBin[i, j] <- 0)
		}
	}

	#Creating range overlap matrix
	overQuant <- matrix(ncol = length(spp), nrow = length(spp))
	rownames(overQuant) <- colnames(overQuant) <- spp

	setScale(1e+09) #era 8
	for (i in 1:length(spp)) {
		for (j in 1:length(spp)) {
			sh_i <- readShapePoly(spp[i])
			sh_j <- readShapePoly(spp[j])
			ifelse(gIntersects(sh_i, sh_j), overQuant[i, j] <- gArea(gIntersection(sh_i, sh_j)), overQuant[i, j] <- 0)
		}
	}
	
	res<-list()
	res$spp<-spp
	res$overBin<-overBin
	res$overQuanto<-overQuant
	res$rangeSizes<-diag(overQuant)
	res$bt <- sort(betweenness(overBin, gmode = "graph"), decreasing = TRUE)
	res$cl <- sort(closeness(overBin, gmode = "graph"), decreasing = TRUE)
	res$dg <- sort(degree(overBin, gmode = "graph"), decreasing = TRUE)
	res
}

