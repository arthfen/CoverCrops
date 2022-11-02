######### Supplementary Material - From regional to parcel scale: a high-resolution map of cover crops across Europe combining satellite data with statistical surveys
####### This file contains the method for parameter estimation
####### Author: Arthur Nicolaus Fendrich
#########

library(terra)
library(mgcv)
library(snow)
library(Rcpp)
library(numDeriv)

# Some preliminary variables and functions for the method

np <- 5 # Number of model parameters
nodes <- 20 # Number of nodes for parallel processing
nr <- nrow(rast('100m/europe.tif')) # Number of rows in the raster file
blocks <- 3000 # Number of blocks for parallel processing

# Loads auxiliary functions written in C++
pw <- getwd()
sourceCpp(paste0(pw, '/utils/matsums.cpp'), cacheDir = paste0(pw, '/tmp'))
source('utils/PredictMat_new.R')

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
cellFromRow <- function(object, rownr) {
	cols <- rep(1:ncol(object), times=length(rownr))
	rows <- rep(rownr, each=ncol(object))	
	cellFromRowCol(object, rows, cols)
}
coordinates <- function(object) {
	xyFromCell(object, cell=1:ncell(object))
}

# Loads the observed variable Y
df.Y <- read.csv('data_percent.csv')
area <- read.csv('data_arablearea.csv')
df.Y$cover_crop_area <- df.Y$cover_crop_percent/100 * area[match(df.Y$id, area$nuts2), 'arable_area_ha']
df.Y <- df.Y[,c('id', 'cover_crop_area')]
df.Y[is.na(df.Y)] <- 0

# Defines the observed variable Y
g <- function(x) 1/(1 + exp(x)) # The link function varies from 0 (no cover crop) to 1 (1 ha - i.e., the pixel size - of cover crop)
dg <- function(x) {
	v = -exp(x)/(exp(x) + 1)^2
	v[is.na(v)] <- 0
	v
}

# From here on, we set-up the environment by creating reference smoothers
# Alternatively, you can skip lines 62-120 and directly load the smoother used in the work in Line 121
for(j in 1:31) assign(paste0('n', j), rast(paste0('100m/IndexesSAll-', j, '.tif')))
block <- rast('100m/europe.tif')
uso <- rast('100m/Corine2018_warp100m-Arable3.tif')

# The areas and land cover
tmpr <- c(block, uso)
tmpr <- spatSample(tmpr, 1e4, method = 'regular', as.raster = TRUE)
tmp <- data.frame(coordinates(tmpr[[1]]), readValues(tmpr, mat = TRUE))
colnames(tmp) <- c('x1', 'x2', 'block', 'uso')

# The index
ndf <- c(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31)
ndf <- spatSample(ndf, 1e4, method = 'regular', as.raster = TRUE)
ndf <- data.frame(readValues(ndf, mat = TRUE)/10000)
colnames(ndf) <- 1:31
ndf[ndf >= 0] <- NA
ndf <- 10^(ndf)
ndf[is.na(ndf)] <- 0

# The coordinates
nco1 <- data.frame(x1 = tmp$x1)
for(j in 2:31) nco1[,paste0('x', j)] <- nco1[,1]
nco2 <- data.frame(x1 = tmp$x2)
for(j in 2:31) nco2[,paste0('x', j)] <- nco2[,1]
colnames(nco1) <- colnames(nco2) <- 1:31
tmp$x1 <- tmp$x2 <- NULL

# The months
mo <- data.frame(rep(1, nrow(ndf)), 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 ,25, 26, 27, 28, 29, 30, 31)
colnames(mo) <- 1:31

# Merge all and remove non-arable land
tmp$cr <- as.matrix(ndf)
tmp$month <- as.matrix(mo)
tmp$x1 <- as.matrix(nco1)
tmp$x2 <- as.matrix(nco2)
tmp$block <- ifelse(tmp$block == 0, NA, tmp$block)
tmp$L <- as.matrix(1 + 0 * ndf)
tmp$L[ndf == 0] <- 0
tmp <- na.omit(tmp)
rm(tmpr, mo, nco1, nco2, ndf)

# Adds minmax values
minmax <- readRDS('output/0minmax.rds') # This is a pre-calculed file with the range of each variable
na <- which.min(tmp$x1); tmp[na, 'x1'] <- minmax[[1]][1,1]; tmp[na, 'x2'] <- minmax[[1]][1,2]
na <- which.max(tmp$x1); tmp[na, 'x1'] <- minmax[[1]][2,1]; tmp[na, 'x2'] <- minmax[[1]][2,2]
na <- which.min(tmp$x2); tmp[na, 'x1'] <- minmax[[2]][1,1]; tmp[na, 'x2'] <- minmax[[2]][1,2]
na <- which.max(tmp$x2); tmp[na, 'x1'] <- minmax[[2]][2,1]; tmp[na, 'x2'] <- minmax[[2]][2,2]
mi <- min(sapply(1:31, function(i) 10^(minmax[[i+2]][1]/10000)))
ma <- 0.5 # Here we use the 0.99 quantile of the CR value instead of the maximum, because the distribution is very skewed
for(i in 1:31){
        tmp$cr[which(tmp$cr[,i] < mi), i] <- mi
#	tmp$cr[which.max(tmp$cr[,i]), i] <- 10^min(0, minmax[[i+2]][2]/10000)
        tmp$cr[which(tmp$cr[,i] > ma), i] <- ma
}
rm(minmax)

# NOTE: for the linear functional term 'month', missing values are naturally treated as zero because this is the NoData value on the 'cr' raster file
sm1 <- smoothCon(te(x1, x2, month, cr, by = L, bs = c('ps', 'ps', 'ps', 'ps'), m = 0, k = c(9, 8, 8, 8)), data = tmp, absorb.cons = TRUE); sm1[[1]]$X <- NULL
#sm1 <- readRDS('output/1arable_sm1.rds')
bd <- 1 + sm1[[1]]$df
S1 <- S2 <- S3 <- S4 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]
S3[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[3]]
S4[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[4]]
#saveRDS(sm1, 'output/1arable_sm1.rds') # Be careful not to overwrite the original file
sm1[[1]]$S <- NULL
for(j in 1:length(sm1[[1]]$margin)) sm1[[1]]$margin[[j]]$S <- sm1[[1]]$margin[[j]]$X <- NULL
rm(tmp)

# Splits the raster in chunks and randomizes to avoid overload
bss <- bssize(nr, blocks)
set.seed(1)
od <- sample(1:bss$n)
bss$row <- bss$row[od]
bss$nrows <- bss$nrows[od]
print(paste0('bss: ', bss$n))
d <- bss$n

# Now, runs the method
# Declares a function to generate the aggregated matrices
ag.m <- function(i, b){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

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
	tmp[tmp[,'block'] == 0, 'block'] <- NA

	na <- getRowIndex_NA(tmp)
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL, NULL))

	tmp <- cbind(xyFromCell(uso, cellFromRow(uso, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1))), tmp)
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
	ndf[ndf < mi] <- 0
	ndf[ndf > ma] <- 0
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
	if(length(na) > 0){
		tmp <- list(L = L[-na,,drop=FALSE], cr = ndf[-na,,drop=FALSE], month = mo[-na,,drop=FALSE], x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na])
	} else {
		tmp <- list(L = L, cr = ndf, month = mo, x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'])
	}
	rm(L, ndf, mo)

	# The ordering of groups in the aggregation function
	groups <- sort(unique(tmp$block))

	# Generates the aggregated matrices using auxiliary C++ code
	tmpgam <- FastPredictMat2_agg(sm1[[1]], s.margin = c(1, 2), t.margin = c(3, 4), tmp, b = b[[1]], qrQ = qrQ)
	rm(tmp)

	# This is equivalent to: return(list(tmp_Agz., tmp_AGX, tmp_MMt, groups))
	print(paste0('sent, i = ', i, ' : ', Sys.time()))

	return(list(tmpgam[['gz.']], tmpgam[['AGX']], tmpgam[['MMt']], groups))
}

# Creates the cluster
cl <- makeCluster(nodes)
clusterExport(cl, list('bss', 'sm1', 'df.Y', 'ag.m', 'g', 'dg', 'mi', 'ma', 'pw', 'cellFromRow', 'coordinates', 'FastPredictMat2_agg', 'S1', 'S2', 'S3', 'S4', 'np'))
clusterApply(cl, seq_along(cl), function(i) sink(paste0('output/worker_', i)))
clusterEvalQ(cl, {
	library(terra)
	library(mgcv)
	library(Rcpp)
	library(numDeriv)
	sourceCpp(paste0(pw, '/utils/matsums.cpp'), cacheDir = paste0(pw, '/tmp/'), verbose = TRUE)

	qrQ <- matrix(1)
	if(!is.null(attr(sm1[[1]], 'qrc'))){
		if(attr(sm1[[1]], 'nCons') > 0){
			qrQ <- t(qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE))
		}
	}

	return(TRUE)
})

# While we don't calculate the analytical gradient explicitly, we use the numerical one
grfitCpp <- function(x, ...) grad(fitCpp, x, ..., method = 'simple')
clusterExport(cl, list('grfitCpp'))

# Here we define a function to iterate over the whole raster
out.fn <- function(b){
	# Creates the auxiliary matrices and proceed
	b <- list(matrix(b, ncol = 1))

	Agz. <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	Agz.[,1] <- df.Y[,1]

	X_ <- matrix(0, nrow = dim(df.Y)[1], ncol = bd)
	diagMMt <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	diagMMt[,1] <- df.Y[,1]

	# Sends different tasks to each node
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = b), tag = i)
	}
	for (i in 1:d){
		k <- recvOneData(cl)

		ni <- nodes + i
		if (ni <= d) {
			sendCall(cl[[k$node]], ag.m, list(i = ni, b = b), tag = ni)
		}
		k <- k$value$value

		if(i == 1 | i %% round(d/4) == 0) print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

		if(!is.null(k[[4]])){
			groups <- k[[4]]

			v.out <- k[[1]]
			v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
			Agz.[,2] <- Agz.[,2] + v.out

			v.out <- k[[2]]
			v.out <- rbind(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1),]
			X_ <- X_ + as.matrix(v.out)

			v.out <- k[[3]]
			v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
			diagMMt[,2] <- diagMMt[,2] + v.out

			rm(v.out, k)
		}
	}
	print(paste0('received: all - ag.m - ', Sys.time()))

	Y_ <- df.Y[,2] - Agz.[,2]
	b <- b[[1]]

	# Generate the auxiliary matrices for the log-likelihood calculation
	diagL <- 0 * diagMMt[,2]
	ql <- diagMMt[,2] > 0
	diagL[ql] <- sqrt(1/diagMMt[ql,2])
	y. <- diagL * Y_
	X. <- diagL * X_
	tp <- qrCpp(X.)
	Qty. <- t(tp$Q) %*% y.
	R <- tp$R
	rm(tp)
	r <- crossprod(y.) - crossprod(Qty.)
	N <- sum(ql)

	er <- crossprod(diagL * Y_)
	print(er)

	return(list(Qty., R, r, N, b))
}

maxit <- 50 # Maximum number of steps
ml <- rep(NA, maxit)
n.dif <- 1

mtcs <- rep(list(NA), maxit)
it <- 0
set.seed(1)
b0 <- rnorm(bd, 0, 0.1) # Initial guess for the predicted random effects, generated randomly here
st <- rep(0, np)

while(n.dif > 1e-3){
	it <- it + 1

	print(paste0('iteration #', it))

	ll.t1 <- out.fn(b=b0)

	mtcs[[it]] <- ll.t1
	mtcs[[it]][[6]] <- st
	saveRDS(mtcs, 'output/1mtcs.rds')
	if(it > 1) saveRDS(mtcs[[it-1]], 'output/1mtcs_sol.rds')

	print(paste0('starting optimization: ', Sys.time()))
	
	clusterExport(cl, list('ll.t1'))
	l <- parLapply(cl, 1:10, function(i){
                set.seed(i)
                optim(rnorm(np), fn=fitCpp, gr=grfitCpp, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], S1=S1, S2=S2, S3=S3, S4=S4, method='BFGS')
        })
	vs <- sapply(l, function(i) i$value)
	print(vs)
	l <- l[[which.min(vs)]]

	ml[it] <- as.numeric(l$value)
	if(it > 2){
		if(ml[it] > ml[it-1]) break
	}

	bnew <- hbCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], S1=S1, S2=S2, S3=S3, S4=S4, w=1)
	gg <- crossprod(grfitCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], S1=S1, S2=S2, S3=S3, S4=S4))

	n.dif <- crossprod(b0 - bnew)
	st <- l$par
	b0 <- bnew

	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('norm gr @ minimum: ', gg))
	print(paste0('ll: ', round(ml[it], 4)))
	print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
	print(paste0('s2: ', round(exp(st[1]), 8)))
}
stopCluster(cl)
quit('no')

