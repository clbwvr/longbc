#' Biclustering Multivariate Longitudinal Data
#' 
#' \code{longbc} estimates the underlying biclustering structure and coefficients of an input data set. This function calls \code{cppadmm} for performing most of the computation.
#' longbc = function(dat, lambdas, V=NULL, tol.abs=1e-3, tol.rel = 1e-3, niter = 20, rho = 1, tausd = 1/4, XtXrhoDeltai=NULL, verbose=TRUE){
#' @param dat The input data.frame with columns "subj","feat","t", and "y" for the subject, feature, time, and response variables, respectively.
#' @param lambdas A vector of regularization parameter values
#' @param V The plug-in estimator for the covariance matrix of the random intercepts. If NULL, a compound symmetric correlation structure will be estimated from the data.
#' @param tol.abs The absolute convergence tolerance
#' @param tol.rel The relative convergence tolerance
#' @param niter The maximum number of iterations
#' @param rho ADMM parameter
#' @param tausd Thresholding parameter
#' @param XtXrhoDeltai Precomputed sparse matrix
#' @param verbose Flag for logging
#' @param trace Flag for objective function tracking
#' @export
#' @author Caleb Weaver
longbc = function(dat, lambdas, V=NULL, tol.abs=1e-2, tol.rel = 1e-3, niter = 50, rho = 1, tausd = 1/4,XtXrhoDeltai=NULL,verbose=TRUE,trace=TRUE){
  
  # Data validation
  if(!all(c("subj","feat","t","y") %in% colnames(dat))){
    stop("dat must contain columns subj, feat, t, and y.")
  }
  
  # Translate responses
  agg = aggregate(data=dat[dat$t==0,], y~feat, mean)
  colnames(agg)[2] = "featmean"
  dat = merge(dat, agg, by="feat")
  dat$y = dat$y-dat$featmean
  dat = dat[,!(colnames(dat) %in% c("featmean","featsd"))]
  dat = dat[order(dat$subj, dat$feat, dat$t),]
  
  # Get input data  
  subj = dat$subj
  feat = dat$feat
  t = dat$t
  y = dat$y
  N = nrow(dat)
  n = length(unique(subj))
  p = length(unique(feat))
  m = length(unique(t))
  g = expand.grid("subj"=unique(subj), "feat"=unique(feat))
  g = g[order(g$subj, g$feat),]
  idx = factor(paste0(g$subj,"_",g$feat))
  unit = factor(paste0(subj,"_",feat))
  unit = factor(unit,levels=unique(unit))
  dat$unit = unit
  X = model.matrix(~unit:t-1)
  colnames(X) = levels(unit)
  X = Matrix(X,sparse=TRUE) 
  b0 = mean(y)
  
  # Initial parameter estimate
  XtX = Matrix::crossprod(X,X)
  Xty = Matrix::crossprod(X,y)
  ols = as.numeric(Matrix::solve(XtX, Xty))
  Xty=as.numeric(Xty)
  dat$pred = as.numeric(X%*%ols)
  dat$res = dat$y - dat$pred
  names(ols) = paste0("beta_",idx)
  gols = expand.grid(unique(subj),unique(feat))
  gols = gols[order(gols[,1]),]
  gols$ols = ols
  colnames(gols) = c("subj","feat","ols")
  Betahat = as.matrix(dcast(data = gols, subj ~ feat, value.var = "ols")[,-1])
  
  # Get Delta matrix
  D = matrix(0,n*p,n*p)
  A = -diag(p)
  for(i in 1:(n-1)){
    for(ip in (i+1):n){
      D[g$subj==i,g$subj==ip] = A
    }
  }
  diag(D) = n-1
  Dt = matrix(0,n*p,n*p)
  for(i in 1:(n)){
    Dt[g$subj==i, g$subj==i] = -1
  }
  diag(Dt) = p-1
  Delta = D + Dt
  
  # Get XtXrhoDeltai matrix
  if(is.null(XtXrhoDeltai)){
    XtXrhoDelta = XtX + rho * Delta
    XtXrhoDeltai = cppmatinv(as.matrix(XtX + rho * Delta))
    XtXrhoDeltai = Matrix(XtXrhoDeltai,sparse=TRUE)
    XtXrhoDeltai = as(XtXrhoDeltai, "dgCMatrix")
  }
  
  #Get distances
  subjdist = matrix(NA, n, n)
  for (i in 1:(n - 1)) {
    for (ip in (i + 1):n) {
      subjdist[i,ip] = cppnrm2(ols[g$subj==i] - ols[g$subj==ip])
    }
  }
  subjdist = cppuppertosym(subjdist)
  
  featdist = matrix(NA, p, p)
  for (i in 1:(p - 1)) {
    for (ip in (i + 1):p) {
      featdist[i,ip] = cppnrm2(ols[g$feat==i] - ols[g$feat==ip])
    }
  }
  featdist = cppuppertosym(featdist)
  
  # Get weights
  wsubj = 1/subjdist
  wfeat = 1/featdist
  wsubj = wsubj * sqrt(n) /  sum(wsubj,na.rm = TRUE)
  wfeat = wfeat * sqrt(p) /  sum(wfeat,na.rm = TRUE)
  
  # Get weight lists
  wsubjlist = list()
  for(i in 1:(n-1)){
    wsubjlist[[i]] = list()
    for(j in (i+1):n){
      wsubjlist[[i]][[j]] = wsubj[i,j]
    }
  }
  wfeatlist = list()
  for(i in 1:(p-1)){
    wfeatlist[[i]] = list()
    for(j in (i+1):p){
      wfeatlist[[i]][[j]] = wfeat[i,j]
    }
  }
  
  # Iterations  
  a.init = ols
  if(verbose) print("Beginning iterations.")
  
  # Lambda loop
  nlambda = length(lambdas)
  fits = list()
  for (li in 1:nlambda) {
    lambda = lambdas[li]
    if(verbose) print(paste0("   lambda = ",round(lambda,3)))
    wsubj = unlist(wsubjlist)
    wfeat = unlist(wfeatlist)
    
    wlrsubj = lapply(wsubjlist, function(u) lapply(u, function(v) v*lambda/rho))
    wlrfeat = lapply(wfeatlist, function(u) lapply(u, function(v) v*lambda/rho))
    
    # Initialize estimates
    if(li==1){
      a = a.init 
    }
    blist = cppsubju(a, g$subj, n)
    clist = cppfeatu(a, g$feat, p)
    deltalist = lapply(blist, function(u) lapply(u,function(v) rep(0, length(v))))
    etalist = lapply(clist, function(u) lapply(u,function(v) rep(0, length(v))))
    u1 = unlist(blist)
    u2 = unlist(clist)
    u = c(u1,u2)
    
    # ADMM Updates
    l = cppadmm(a,blist,clist,deltalist,etalist,g$subj,g$feat,n,p,X,y,lambda,wsubj,wfeat,wlrsubj,wlrfeat,Xty, XtXrhoDeltai, rho, niter, tol.rel, tol.abs, trace, verbose)
    
    a = l[[1]]
    losses = l[[2]]
    names(a) = idx
    g$beta = a
    Betahat = dcast(data = g, subj ~ feat, value.var = "beta")
    Betahat = as.matrix(Betahat[, 2:ncol(Betahat)])
    
    # Thresholding
    vlist =  list()
    taumat = matrix(NA,n,n)
    for(i in 1:(n-1)){
      vlist[[i]] = list()
      for(j in (i+1):n){
        vlist[[i]][[j]] = Betahat[i,] - Betahat[j,]
        taumat[i,j] = cppnrm2(vlist[[i]][[j]])
      }
    }
    tau = sd(c(taumat),na.rm=TRUE)
    
    vtlist = list()
    tautmat = matrix(NA,p,p)
    for(i in 1:(p-1)){
      vtlist[[i]] = list()
      for(j in (i+1):p){
        vtlist[[i]][[j]] = Betahat[,i] - Betahat[,j]
        tautmat[i,j] = cppnrm2(vtlist[[i]][[j]])
      }
    } 
    taut = sd(c(tautmat),na.rm=TRUE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if(taumat[i,j] < tau * tausd) taumat[[i,j]]=0
      }
    }
    
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(tautmat[i,j] < taut * tausd) tautmat[i,j]=0
      }
    }
    
    # Bicluster assignment
    hcs = hclust(as.dist(t(taumat)))
    cs = cutree(hcs,h = 0)
    dfcs = data.frame(subj=unique(subj), cs=cs)
    
    hcf = hclust(as.dist(t(tautmat)))
    cf = cutree(hcf,h = 0)
    dfcf = data.frame(feat=unique(feat), cf=cf)
    df = merge(g,dfcs,by="subj")
    df = merge(df,dfcf,by="feat")
    df = df[order(df$subj, df$feat),]
    cb = as.numeric(factor(paste(df$cs,df$cf))) 
    df$cb = cb 
    
    d1 = dat
    d1 = merge(d1[,c("subj","feat","t","y")],df[,c("subj","feat","cb")], by=c("subj", "feat"))
    res = lm(y ~ t * cb - 1,data=d1)
    bic = BIC(res)
    Muhat = 0 * Betahat
    for(i in unique(cs)){
      for(j in unique(cf)){
        Muhat[cs==i, cf==j] = mean(Betahat[cs==i, cf==j])
      }
    }
    
    fits[[li]] = list(
      df=df, 
      Betahat = Betahat, 
      Muhat = Muhat,
      losses = losses,
      cs = cs,
      cf = cf,
      cb = cb,
      lambda = lambda,
      bic = bic
    )
    
  }
  
  obj = list()
  obj[['fits']] = fits
  obj[["XtXrhoDeltai"]] = XtXrhoDeltai
  bics = sapply(fits, function(u) u$bic)
  obj[['dat']] = dat
  sd = diff(range(bics[-1]))/5
  obj[['bestfit']] = max(which(bics == min(bics)))
  obj[["bics"]] = bics
  class(obj) = "longbc"
  
  return(obj)
}

#' @export
ari.longbc = function(obj,truth){
  cblist = lapply(obj$fits, function(u) u$df[,c("subj","feat","cb")])
  aris = sapply(cblist, function(u){
    truth = truth[,c("subj","feat","cb")]
    df1 = merge(u, truth, by=c("subj","feat"))
    adjustedRandIndex(df1$cb.x,df1$cb.y)
  })
  return(aris)
}

#' @export
bic.longbc = function(obj){
  cblist = lapply(obj$fits, function(u) u$df[,c("subj","feat","cb")])
  dat1 = obj$dat
  bics = sapply(cblist, function(u){
    dat2 = dat1[,c("subj","feat","t","y")]
    dat2$unit = factor(paste0(dat2$subj,"_",dat2$feat))
    dat3 = merge(dat2, u, by=c("subj","feat"))
    res = lm(y ~ (t|cb),  data=dat3)
    BIC(res)
  })
  return(bics)
}