######### Supplementary Material - From regional to parcel scale: a high-resolution map of cover crops across Europe combining satellite data with statistical surveys
####### This file must be run to convert the exported RDS files to a single compressed TIFF
####### Author: Arthur Nicolaus Fendrich
#########

library(raster)

nr <- nrow(raster('../100m/europe.tif'))
blocks <- 2500
bssize <- function(nr, nb) {
	size <- ceiling(nr/nb)
	nb <- ceiling(nr / size)

	row <- (0:(nb-1))*size + 1
	nrows <- rep(size, length(row))
	dif <- nb * size - nr 
	nrows[length(nrows)] = nrows[length(nrows)] - dif
	return(list(row=row, nrows=nrows, n=nb))
}

bss <- bssize(nr, blocks)

n1 <- raster('../100m/IndexesSAll-1.tif')
out_mean <- raster(n1)
out_mean <- raster::writeStart(out_mean, filename = '3output_mean', datatype = 'INT4U', format='GTiff', options = "COMPRESS=LZW")
for (i in 1:bss$n) {
	k <- readRDS(paste0('tifs/v_', i))
	print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

	out_mean <- raster::writeValues(out_mean, k, bss$row[i])
}
out_mean <- raster::writeStop(out_mean)

