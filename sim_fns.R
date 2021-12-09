source("bicluster_fns.R")
library(sparseBC)
library(cvxbiclustr)
library(Rcpp)
library(dplyr)
library(reshape2)
library(stats)
library(mclust)
library(unikn)
library(gplots)

#' @export
parmf = function(ns = c(20,100), ps = c(10,30), snrs = c(.5, 1),ms=c(5), nbs = c(9, 16),grids=c(TRUE,FALSE),exchs=c(TRUE,FALSE)){
  parammat = expand.grid(n=ns, p=ps, snr=snrs, nb=nbs,m=ms, grid=grids,exch=exchs)
  parammat = parammat[order(parammat$n,parammat$p,parammat$snr,parammat$nb,parammat$grid,parammat$exch),]
  paramlist = lapply(1:nrow(parammat), function(i) parammat[i,])
  return(paramlist)
}

### Simulation function
#' @export
datf = function(parm){
  n = parm$n
  p = parm$p
  snr = parm$snr
  nb = parm$nb
  m = parm$m
  exch = parm$exch
  grid = parm$grid
  salpha = 0
  if(exch){ #exchangable
    R = matrix(.5, p, p)
    diag(R) = 1
    P = diag(sqrt(salpha),p,p) %*% R %*% diag(sqrt(salpha),p,p) 
    V = kronecker(diag(n),P)
    V = Matrix(V, sparse=TRUE)
    alpha = Rfast::rmvnorm(100, rep(0,n*p), V)
  } else{ # independent
    # R = matrix(0, p, p)
    # diag(R) = 1
    # P = diag(sqrt(salpha),p,p) %*% R %*% diag(sqrt(salpha),p,p) 
    # V = kronecker(diag(n),P)
    # V = Matrix(V, sparse=TRUE)
    # alpha = c(Rfast::rmvnorm(1, rep(0,n*p), V))
    alpha = rep(0,n*p)
  }
  beta = sample((-10):10,nb,FALSE)
  subj = 1:n
  feat = 1:p
  dfs = data.frame(subj=subj, cs=sample(1:sqrt(nb),n,TRUE))
  dff = data.frame(feat=feat, cf=sample(1:sqrt(nb),p,TRUE))
  g = expand.grid(subj=subj, feat=feat)
  g = merge(g, dfs, by="subj")
  g = merge(g, dff, by="feat")
  g = g[order(g$subj, g$feat),]
  g$unit = factor(paste0(g$subj,"_",g$feat))
  g$cb = as.numeric(as.factor(paste(g$cs, g$cf)))
  CB = dcast(data=g, subj ~ feat, value.var="cb")
  CB = as.matrix(CB[,-1])
  alphaij = alpha[as.numeric(g$unit)]
  betaij = beta[g$cb]
  g$alpha = alphaij
  g$beta = betaij
  Beta = dcast(data=g, subj ~ feat, value.var="beta")
  Beta = as.matrix(Beta[,-1])
  
  # Create design matrix
  if(grid) {
    dat = expand.grid(subj=1:n, feat=1:p, t=0:(m-1)); 
    dat = dat[order(dat$subj, dat$feat, dat$t),]
  } else{
    subjs = feats = ts = c()
    for(i in 1:n){
      mi = sample(c(3,4,5,6),1,prob=c(1,.9,.8,.7))
      ti = 0:(mi-1)
      for(j in 1:p){
        subjs = c(subjs, rep(i, mi))
        feats = c(feats, rep(j, mi))
        ts = c(ts, ti)
      }
    }
    dat = data.frame(subj = subjs, feat=feats, t=ts)
  }
  
  datx = dat
  datx$row = 1:nrow(datx)
  datx$int = 1
  X = dcast(data=datx, row ~ subj+feat, value.var="t",fill = 0)[,-1]
  X = Matrix(as.matrix(X), sparse=TRUE)
  Zi = kronecker(diag(p), rep(1,m))
  Z = dcast(data=datx, row ~ subj+feat, value.var="int",fill = 0)[,-1]
  Z = Matrix(as.matrix(Z), sparse=TRUE)
  ey = as.numeric(Z%*%alphaij + X%*%betaij)
  N = length(ey)
  white = rnorm(N)
  signal = var(ey)
  k = sqrt(signal) / snr * var(white)
  e = k * white
  y = ey + e
  
  # Get Delta matrix
  D = matrix(0,n*p,n*p)
  A = -diag(p)
  for(i in 1:(n-1)){
    for(ip in (i+1):n){
      D[g$subj==i,g$subj==ip] = A
    }
  }
  D = upper.to.sym(D)
  diag(D) = n-1

  Dt = matrix(0,n*p,n*p)
  for(i in 1:(n)){
    Dt[g$subj==i, g$subj==i] = -1
  }
  diag(Dt) = p-1

  Delta = D + Dt
  
  rho = 1
  XtXrhoDelta = Matrix::crossprod(X,X) + rho * Delta
  XtXrhoDeltai = matrix(solve(XtXrhoDelta),n*p,n*p)
  XtXrhoDeltai[matrix(XtXrhoDelta==0,n*p,n*p)] = 0
  XtXrhoDeltai = Matrix(XtXrhoDeltai,sparse=TRUE)
  XtXrhoDeltai = as(XtXrhoDeltai, "dgCMatrix")
  dat$ey = ey
  dat$y = y
  agg = aggregate(data=dat[dat$t==0,], y~feat, mean)
  colnames(agg)[2] = "featmean"
  dat = merge(dat, agg, by="feat")
  dat$y = dat$y-dat$featmean
  
  agg = aggregate(data=dat[dat$t==0,], y~feat, sd)
  colnames(agg)[2] = "featsd"
  dat = merge(dat, agg, by="feat")
  dat$y = dat$y/dat$featsd
  
  dat = dat[,!(colnames(dat) %in% c("featmean","featsd"))]
  dat = merge(dat, g, by=c("subj","feat"))
  dat = dat[order(dat$subj, dat$feat, dat$t),]
  
  # # Initial parameter estimate
  XtX = Matrix::crossprod(X,X)
  Xty = Matrix::crossprod(X,y)
  ols = as.numeric(Matrix::solve(XtX, Xty))
  g = dat[,c("subj","feat","cs","cf","cb")]
  g = g[!duplicated(g),]
  g$ols = ols
  
  return(list(parm=parm, dat=dat, XtXrhoDeltai=XtXrhoDeltai,V = NULL,ols= NULL))
}

dofnaive = function(datres){
  parm = datres[[1]]
  dat = datres[[2]]
  XtXrhoDeltai = datres[[3]]
  V = datres[[4]]
  ols = datres[[5]]
  betaij = dat$beta
  g = dat[,c("subj","feat","cs","cf","cb","beta")]
  g = g[!duplicated(g),]
  g$ols = ols
  betaij = g$beta
  
  n = parm$n
  p = parm$p
  snr = parm$snr
  nb = parm$nb
  snb = sqrt(nb)
  exch = parm$exch
  grid = parm$grid
  
  ### 
  # OLS 
  Beta = as.matrix(dcast(data=g, subj ~ feat, value.var="beta")[,-1])
  Betahat = as.matrix(dcast(data=g, subj ~ feat, value.var="ols")[,-1])
  ari.ols = 0
  mse.ols1 = mean((ols - betaij)^2)
  mse.ols = mean((Beta - Betahat)^2)
  print(mse.ols)
  print(mse.ols1)
  # matshow(Betahat)
  nb.ols = nrow(g)
  ###
  
  ### 
  # DCT 
  # distr = dist(Betahat)
  # distc = dist(t(Betahat))
  # hr = hclust(distr)
  # hc = hclust(distc)
  # r = dynamicTreeCut::cutreeHybrid(hr, distM = as.matrix(distr),verbose=0)$labels
  # c = dynamicTreeCut::cutreeHybrid(hc, distM = as.matrix(distc),verbose=0)$labels
  # g1=expand.grid(c,r)
  # g1$cb = as.numeric(factor(paste(g1$Var1,g1$Var2)))
  # ari0 = adjustedRandIndex(g1$cb, g$cb)
  # nb0 = length(unique(g1$cb))
  #matshow(Betahat[order(r),order(c)])
  ###
  
  ### 
  # Sparse BC
  set.seed(NULL)
  res1 = spbc.cv(Betahat, 1:min(n-1, 2*snb),1:min(p-1, 2*snb),0)
  mse.bc = mean((res1[[2]] - Beta)^2)
  ari1 = adjustedRandIndex(res1[[1]]$cb, g$cb)
  nb1 = length(unique(res1[[1]]$cb))
  gcs = res1[[1]][,c("subj","cs")]
  cs1 = gcs[!duplicated(gcs),"cs"]
  gcf = g[,c("feat","cf")]
  cf1 = gcf[!duplicated(gcf),"cf"]
  
  ### 
  # COBRA
  # A = Betahat
  # A = A - sum(A)/length(A)
  # A = A/norm(A,"f")
  # wts = cvxbiclustr::gkn_weights(A)
  # w_row = wts$w_row
  # w_col = wts$w_col
  # E_row = wts$E_row
  # E_col = wts$E_col
  # #### Initialize path parameters and structures
  # nGamma = 10
  # ## Generate solution path
  # sol = NULL
  # mg = 0
  # gammaSeq <- 10**seq(.5,1,length.out=nGamma)
  # sol <- hush({cobra_validate(A,E_row,E_col,w_row,w_col,gammaSeq)})
  # sink()
  # verr = sol$validation_error
  # ix = max(which(verr == min(verr)))
  # U = sol$U[[ix]]
  # r = sol$groups_row[[ix]]$cluster
  # c = sol$groups_col[[ix]]$cluster
  # g1 = expand.grid(c,r)
  # cb = as.numeric(factor(paste(g1[,1], g1[,2])))
  # ari2 = adjustedRandIndex(cb, g$cb)
  # nb2 = length(unique(cb))
  return(list(dat=dat, nb=nb1, ari=ari1))
}

dofprop = function(datres,maxlam){
  parm = datres[[1]]
  dat = datres[[2]]
  XtXrhoDeltai = datres[[3]]
  V = datres[[4]]
  ols = datres[[5]]
  betaij = dat$beta
  g = dat[,c("subj","feat","cs","cf","cb","beta")]
  g = g[!duplicated(g),]
  Beta = dcast(data=g, subj ~ feat, value.var="beta")
  Beta = as.matrix(Beta[,-1])
  g$ols = ols
  betaij = g$beta
  
  n = parm$n
  p = parm$p
  snr = parm$snr
  nb = parm$nb
  snb = sqrt(nb)
  exch = parm$exch
  grid = parm$grid
  
  # Longitudinal BC
  lambdas = seq(0,maxlam,length.out=15)
  fit = longbc:::longbc(dat = dat, lambdas = lambdas, tausd=1/4, rho=1, tol.abs=1e-2, tol.rel=1e-2,niter = 20, XtXrhoDeltai = XtXrhoDeltai,trace=TRUE,nnk=4,wtype="l2",sparsew=FALSE,loud=FALSE)
  aris = longbc:::ari.longbc(fit, g)
  bics = longbc:::bic.longbc(fit)
  best = which.max(aris)
  ari = max(aris)
  nb = length(unique(fit$fits[[best]]$cb))
  return(mean(Beta - fit$fits[[best]]$Betahat)^2)
#  return(list(fit=fit, dat=dat, nb=nb, ari=ari, arifit = aris[which.min(bics)],Betahat = fit[[1]][[best]]$Betahat))
}

#' @export
plotdat = function(dat){
  par(mfrow=c(1,2))
  nb = length(unique(dat$cb))
  datmat = dcast(data=dat, cb + subj + feat  ~ t, value.var="ey")
  matplot(t(datmat[,-c(1,2,3)]), type="l",  col=rainbow(nb)[datmat$cb], lty=1)
  datmat = dcast(data=dat, cb + subj + feat  ~ t, value.var="y")
  matplot(t(datmat[,-c(1,2,3)]), type="l", col=rainbow(nb)[datmat$cb], lty=1)
  par(mfrow=c(1,1))
}

