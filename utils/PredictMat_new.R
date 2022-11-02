######### Supplementary Material - From regional to parcel scale: a high-resolution map of cover crops across Europe combining satellite data with statistical surveys
####### This file contains the row-wise Kronecker product and the corresponding aggregation into matrices, if convenient
####### Author: Zheyuan Li and Arthur Nicolaus Fendrich
#########

## `sm` is the result of `sm <- SmoothCon(...)[[1]]`
## `s.margin` is the ID for spatial margins
## `t.margin` is the ID for temporal margins
## `newdata` is a list (without NA) giving essential covariate values
FastPredictMat2_agg <- function (sm, s.margin, t.margin, newdata, b, qrQ) {
  ## extract margins
  sm.s <- sm$margin[s.margin]
  sm.t <- sm$margin[t.margin]
  by.var <- sm$by
  nb <- dim(newdata[[by.var]])[1]
  ## number of spatial and temporal margins
  n.s <- length(sm.s)
  n.t <- length(sm.t)
  ## total number of margins
  n <- n.s + n.t
  ## deal with spatial margins
  Xs <- vector("list", n.s)
  i <- 1
  while (i <= n.s) {
    Xs[[i]] <- Predict.matrix2(sm.s[[i]], newdata)
    i <- i + 1
  }
  Xs <- tensor.prod.model.matrix(Xs)
  ## then, deal with temporal margins
  nb <- dim(newdata[[by.var]])[2]
  for(j in 1:nb){
    Xt <- vector("list", n.t)
    i <- 1
    while (i <= n.t) {
      Xt[[i]] <- list(newdata[[sm.t[[i]]$term]][,j])
      names(Xt[[i]]) <- sm.t[[i]]$term
      Xt[[i]] <- Predict.matrix2(sm.t[[i]], Xt[[i]])
      i <- i + 1
    }
    Xt <- tensor.prod.model.matrix(Xt)
    if(j == 1) X <- newdata[[by.var]][,j] * Xt
    else X <- X + newdata[[by.var]][,j] * Xt

    rm(Xt)
  }
  ## build model matrix
  ## TODO: very important! now we can not handle centering constraints anymore
  block <- newdata[['block']]

  X <- tpmm_agg_noblock(Xs = Xs, Xt = X, groups = sort(unique(block)), block = block, b = b, qrQ = qrQ)

  ## return X
  X
}

FastPredictMat2_vec_q <- function (sm, s.margin, t.margin, newdata, b, qrQ, type) {
  ## extract margins
  sm.s <- sm$margin[s.margin]
  sm.t <- sm$margin[t.margin]
  by.var <- sm$by
  nb <- dim(newdata[[by.var]])[1]
  ## number of spatial and temporal margins
  n.s <- length(sm.s)
  n.t <- length(sm.t)
  ## total number of margins
  n <- n.s + n.t
  ## deal with spatial margins
  Xs <- vector("list", n.s)
  i <- 1
  while (i <= n.s) {
    Xs[[i]] <- Predict.matrix2(sm.s[[i]], newdata)
    i <- i + 1
  }
  Xs <- tensor.prod.model.matrix(Xs)
  ## then, deal with temporal margins
  nb <- dim(newdata[[by.var]])[2]
  for(j in 1:nb){
    Xt <- vector("list", n.t)
    i <- 1
    while (i <= n.t) {
      Xt[[i]] <- list(newdata[[sm.t[[i]]$term]][,j])
      names(Xt[[i]]) <- sm.t[[i]]$term
      Xt[[i]] <- Predict.matrix2(sm.t[[i]], Xt[[i]])
      i <- i + 1
    }
    Xt <- tensor.prod.model.matrix(Xt)
    if(j == 1) X <- newdata[[by.var]][,j] * Xt
    else X <- X + newdata[[by.var]][,j] * Xt

    rm(Xt)
  }
  ## build model matrix
  block <- newdata[['block']]

  Xv <- tpmm_vec_q_noblock(Xs = Xs, Xt = X, b = b, qrQ = qrQ, type)

  ## return Xv
  Xv
}

