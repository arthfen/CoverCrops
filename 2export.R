######### Supplementary Material - From regional to parcel scale: a high-resolution map of cover crops across Europe combining satellite data with statistical surveys
####### This file contains the codes to export the TIFF files
####### Author: Arthur Nicolaus Fendrich
#########

rm(list=ls())
library(raster)
library(mgcv)
library(snow)
library(parallel)
library(Rmpi)
library(Rcpp)
library(mallinfo)

nodes <- 20
nr <- nrow(raster('100m/europe.tif'))
blocks <- 2500

# Loads auxiliary functions written in C++
pw <- getwd()
sourceCpp(paste0(pw, '/utils/matsums.cpp'), cacheDir = paste0(pw, '/tmp'))
source('utils/PredictMat_new.R')

cl <- makeCluster(nodes, type = 'MPI')
clusterEvalQ(cl, library(terra))
clusterEvalQ(cl, library(mgcv))
clusterEvalQ(cl, library(Rcpp))
clusterEvalQ(cl, library(mallinfo))
clusterExport(cl=cl, list('pw'))
clusterEvalQ(cl, sourceCpp(paste0(pw, '/utils/matsums.cpp'), cacheDir = paste0(pw, '/tmp')))

# Defines auxiliary functions for spatial processing
bssize <- function(nr, nb) {
	size <- ceiling(nr/nb)
	nb <- ceiling(nr / size)

	row <- (0:(nb-1))*size + 1
	nrows <- rep(size, length(row))
	dif <- nb * size - nr 
	nrows[length(nrows)] = nrows[length(nrows)] - dif
	return(list(row=row, nrows=nrows, n=nb))
}
cellFromRow_t <- function(object, rownr) {
	cols <- rep(1:ncol(object), times=length(rownr))
	rows <- rep(rownr, each=ncol(object))	
	cellFromRowCol(object, rows, cols)
}

# Defines the observed variable Y
g <- function(x) 1/(1 + exp(x)) # The link function varies from 0 (no cover crop) to 1 (1 ha - i.e., the pixel size - of cover crop)
dg <- function(x) -exp(x)/(exp(x) + 1)^2

# Loads the reference smoother
sm1 <- readRDS('output/1arable_sm1.rds')
bd <- 1 + sm1[[1]]$df
S1 <- S2 <- S3 <- S4 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]
S3[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[3]]
S4[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[4]]

mtcs = readRDS('output/1mtcs.rds')
it = 4 # The last successful iteration

# The two results
st = mtcs[[it]][[6]]
b0 = mtcs[[it]][[5]]
# The matrices to generate them
Qty=mtcs[[it-1]][[1]]
R=mtcs[[it-1]][[2]]
r=mtcs[[it-1]][[3]]
N=mtcs[[it-1]][[4]]

# Calculates the variance-covariance matrix to generate posterior simulations
rho <- exp(st)
Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4
V <- crossprod(R)/rho[1] + Sl
eigenV <- eigen(V)
sqrtinv <- sqrt(1/eigenV$values)
hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))
hat.beta <- b0
set.seed(1)
Cv <- chol(hat.cov, pivot = TRUE)
Cv <- Cv[,order(attr(Cv, 'pivot'))]
n.rep = 1000 # Number of simulations, can be increased for more precise results
nb <- dim(hat.cov)[1]
br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep)
br <- as.numeric(hat.beta) + br
hat.beta <- br

# Splits the raster in chunks
bss <- bssize(nr, blocks)

# Declares a function to generate the predicted vectors
ag.m <- function(i, b){
	block <- rast('100m/europe.tif')
	uso <- rast('100m/Corine2018_warp100m-Arable3.tif')
	for(j in 1:31) {
		assign(paste0('n', j), rast(paste0('100m/IndexesSAll-', j, '.tif')))
		readStart(get(paste0('n', j)))
	}
	readStart(block)
	readStart(uso)

	## The areas and land cover
	tmp <- cbind(
		readValues(block, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(uso, row = bss$row[i], nrows = bss$nrow[i]))
	colnames(tmp)[1:2] <- c('block', 'uso')
	table(tmp[,'block'])

	# Only France?
#	tmp[!tmp[,'block'] %in% 92:113, 'block'] <- NA

	tmp[tmp[,'block'] == 0, 'block'] <- NA
	na <- getRowIndex_NA(tmp)
	if(length(na) == nrow(tmp)) return(NA * tmp[,1])

	tmp <- cbind(xyFromCell(uso, cellFromRow_t(uso, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1))), tmp)
	colnames(tmp)[1:2] <- c('x1', 'x2')

	## The index
	ndf <- cbind(
		readValues(n1, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n2, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n3, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n4, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n5, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n6, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n7, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n8, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n9, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n10, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n11, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n12, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n13, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n14, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n15, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n16, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n17, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n18, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n19, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n20, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n21, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n22, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n23, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n24, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n25, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n26, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n27, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n28, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n29, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n30, row = bss$row[i], nrows = bss$nrow[i]),
		readValues(n31, row = bss$row[i], nrows = bss$nrow[i]))/10000
	colnames(ndf) <- 1:31
	ndf <- 10^(ndf)
	ndf[ndf >= 0.5] <- 0
	for(j in 1:31) {
		readStop(get(paste0('n', j)))
	}
	readStop(block)
	readStop(uso)

	L <- 1 + 0 * ndf
	L[ndf == 0] <- 0

	## The months
	mo <- matrix(1:31, nrow = nrow(ndf), ncol = 31, byrow = TRUE)

	## Merge all and remove non-arable land
	outval <- NA * ndf[,1]
	if(length(na) > 0){
		tmp <- list(L = L[-na,,drop=FALSE], cr = ndf[-na,,drop=FALSE], month = mo[-na,,drop=FALSE], x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na])
	} else {
		tmp <- list(L = L, cr = ndf, month = mo, x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'])
	}
	rm(L, ndf, mo)

	# Generates the aggregated matrices using auxiliary C++ code
	# Type:
	## 0: median
	## 1: standard deviation
	## 2: 0.05 quantile
	## 3: 0.95 quantile
	tmpgam <- FastPredictMat2_vec_q(sm1[[1]], s.margin = c(1, 2), t.margin = c(3, 4), tmp, b = b[[1]], qrQ = qrQ, type=0)
	outval[-na] <- tmpgam
	rm(tmp)

	return(round(10000 * outval, 0))
}

qrQ <- matrix(1)
if(!is.null(attr(sm1[[1]], 'qrc'))){
	if(attr(sm1[[1]], 'nCons') > 0){
		qrQ <- t(qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE))
	}
}
clusterExport(cl=cl, list('sm1', 'bss', 'cellFromRow_t', 'ag.m', 'g', 'hat.beta', 'qrQ', 'FastPredictMat2_vec_q'))

# Here we define a function to iterate and export the target raster
exporta <- function(){
	l = 0
        for(i in 1:nodes){
                sendCall(cl[[i]], ag.m, list(i = i+l, b = list(hat.beta)), tag = i+l)
        }

        for (i in 1:bss$n) {
                d <- recvOneData(cl)

                if (!d$value$success) {
                        saveRDS(d, 'error.rds')
                        cat('error in number: ', d$value$tag, '\n'); flush.console();
                        stop('cluster error')
                }

                ni <- nodes + i + l
                if (ni <= bss$n) {
                        sendCall(cl[[d$node]], ag.m, list(i = ni, b = list(hat.beta)), tag = ni)
                }

                b <- d$value$tag
		print(paste0('received: ', b))

		saveRDS(d$value$value, paste0('output/tifs/v_', b))
                rm(d)
        }
}
exporta()

# In the end, we have to manually run utils/rds2tif.R to convert the generated RDS files to a single compressed TIF
