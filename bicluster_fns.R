# Helper functions
rc.to.mats = function(X, r, c){
  n = nrow(X)
  p = ncol(X)
  Mu = matrix(NA,n,p)
  for(i in 1:n){
    for(j in 1:p){
      Mu[i,j]=paste0(r[i],"_",c[j])
    }
  }
  Mu=matrix(as.numeric(as.factor(Mu)),n,p)
  K=min(c(max(Mu),50))
  rmat = matrix(FALSE,nrow=n,ncol=K)
  cmat = matrix(FALSE,nrow=K,ncol=p)
  for(k in 1:K){
    for(i in 1:nrow(Mu)){
      rmat[i,k] = as.logical(max(Mu[i,]==k))
    }
    for(j in 1:ncol(Mu)){
      cmat[k,j] = as.logical(max(Mu[,j]==k))
    }
  }
  return(list(rmat=rmat,cmat=cmat))
}

res.to.B = function(X, res){
  n = nrow(X)
  p = ncol(X)
  K = res@Number
  if(K==0) return(matrix(0, n, p))
  if(K==1) cmat=t(res@NumberxCol)
  else cmat = res@NumberxCol
  B = matrix(0,n,p)
  for(i in 1:n){
    for(j in 1:p){
      for(k in 1:K){
        if(res@RowxNumber[i,k] & cmat[k,j]) B[i,j] = k
      }
    }
  }
  return(B)
}

res.to.Mu = function(X, res){
  n = nrow(X)
  p = ncol(X)
  K = res@Number
  if(K==0) return(matrix(0, n, p))
  if(K==1) cmat=t(res@NumberxCol)
  B = res.to.B(X, res)
  Mu = matrix(0,n,p)
  for(k in 1:K){
    Mu[B==k]=mean(X[B==k])
  }
  
  return(Mu)
}

# Jaccard matrix
jmat = function(res1, res2) {
  if(res1@Number==0 | res2@Number==0) mat = matrix(0, ncol = 1, nrow = 1)
  if(res1@Number==1){
    res1@NumberxCol = t(res1@NumberxCol)
  }
  if(res2@Number==1){
    res2@NumberxCol = t(res2@NumberxCol)
  }
  mat = matrix(nrow = res1@Number, ncol = res2@Number)
  rownames(mat) = paste("BC", 1:res1@Number, sep = "")
  colnames(mat) = paste("BC", 1:res2@Number, sep = "")
  for (i in 1:res1@Number) {
    for (j in 1:res2@Number) {
      A = which(res1@RowxNumber[, i] %*% t(res1@NumberxCol[i, ]) > 0)
      B = which(res2@RowxNumber[, j] %*% t(res2@NumberxCol[j, ]) > 0)
      C = intersect(A, B)
      mat[i, j] = length(C)/(length(A) + length(B) - length(C))
    }
  }
  return(mat)
}

# Return recovery and relevance
bc.eval = function(res1, res2, sparse1=TRUE, sparse2=TRUE){
  q1 = ifelse(sparse1, 2,1)
  q2 = ifelse(sparse1, 2,1)
  if(res1@Number==0 | res2@Number==0){
    j1=0
    j2=0
  }
  else{
    mat=jmat(res1, res2)
    j1=mean(apply(mat,2,max))
    j2=mean(apply(mat,1,max))
  }
  return(list(recovery=j1, relevance=j2))
}

# Choose k and r for spBC
sparseBC.choosekr = function (x, k, r, lambda, percent = 0.1, trace = FALSE)  {
  miss = sample(1:(nrow(x) * ncol(x)), nrow(x) * ncol(x), replace = FALSE)
  numberoftimes = 1/percent
  allresults = array(NA, dim = c(numberoftimes, length(k), 
                                 length(r)))
  Cs.init = matrix(NA, nrow = nrow(x), ncol = length(k))
  for (i in 1:length(k)) {
    Cs.init[, i] = kmeans(x, k[i], nstart = 20)$cluster
  }
  Ds.init = matrix(NA, nrow = ncol(x), ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] = kmeans(t(x), r[j], nstart = 20)$cluster
  }
  for (i in 1:numberoftimes) {
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss = x
    missing = miss[1:(nrow(x) * ncol(x))/numberoftimes]
    xmiss[missing] = NA
    xmiss[missing] = mean(xmiss, na.rm = TRUE)
    for (a in 1:length(k)) {
      for (b in 1:length(r)) {
        res = sparseBC(xmiss, k[a], r[b], lambda = lambda, 
                       Cs.init = Cs.init[, a], Ds.init = Ds.init[, 
                                                                 b])$mus
        allresults[i, a, b] = sum((x[missing] - res[missing])^2)
      }
    }
    miss = miss[-1:-(dim(x)[1] * dim(x)[2]/numberoftimes)]
  }
  results.se = apply(allresults, c(2, 3), sd)/sqrt(numberoftimes)
  results.mean = apply(allresults, c(2, 3), mean)
  IndicatorMatrix = 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                              1)] <= results.mean[2:length(k), 2:length(r)] + results.se[2:length(k), 
                                                                                                                         2:length(r)])
  if (max(IndicatorMatrix) == 0) 
    return(list(bestK = max(k), bestR = max(r)))
  RowIndexPlusColIndex = outer(k[-length(k)], r[-length(r)], 
                               "*")
  smallestIndicatorTrue = min(RowIndexPlusColIndex[IndicatorMatrix == 
                                                     TRUE])
  out = which(IndicatorMatrix == TRUE & RowIndexPlusColIndex == 
                smallestIndicatorTrue, arr.ind = TRUE)
  out = matrix(c(k[out[, 1]], r[out[, 2]]), nrow(out), ncol(out))
  temprow = NULL
  tempcol = NULL
  for (i in 1:length(k)) {
    temprow = c(temprow, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempcol = c(tempcol, paste("R = ", r[i], sep = ""))
  }
  rownames(results.se) = temprow
  colnames(results.se) = tempcol
  rownames(results.mean) = temprow
  colnames(results.mean) = tempcol
  return(list(estimated_kr = out, results.se = results.se, 
              results.mean = results.mean))
}

# Function for bicluster analysis
bc.f = function(X, method,alpha=1,gamma=NULL,delta=1,k=1,r=1,k1=1,k2=1,weightsk=5,gamma1=.1,gamma2=.01,...){
  if(method=="ssvd"){
    sink("file.txt") # keep quiet
    res = biclust(X,method=BCssvd(),...)
    sink(); sink();
  }
  else if(method=="plaid"){
    res = biclust(X,method=BCPlaid(), verbose = FALSE,...)
  }
  else if(method=="cc"){
    res = biclust(X,method=BCCC(),number=10,delta=delta,alpha=alpha)
  }
  else if(method=="isa"){
    isa.result = isa(X)
    res = isa.biclust(isa.result)
  }
  else if(method=="dct"){
    distr = dist(X)
    distc = dist(t(X))
    hr = hclust(distr)
    hc = hclust(distc)
    r = cutreeHybrid(hr, distM = as.matrix(distr),verbose=0)$labels
    c = cutreeHybrid(hc, distM = as.matrix(distc),verbose=0)$labels
    mats = rc.to.mats(X,r,c)
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  }  
  else if(method=="spbc"){
    lambda=sparseBC.BIC(X,k=k,r=r,lambda=c(0,10,20,30,40))$lambda
    sol = sparseBC(X,k,r,lambda)
    r = sol$Cs
    c = sol$Ds
    mats = rc.to.mats(X,r,c)
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  }  
  else if(method=="hclust"){
    distr = dist(X)
    distc = dist(t(X))
    hr = hclust(distr)
    hc = hclust(distc)
    r = cutree(hr, k1)
    c = cutree(hc, k2)
    mats = rc.to.mats(X,r,c)
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  }  
  else if(method=="scobra"){
    fit = fit.ssbc(X,gamma1 = gamma1, gamma2=gamma2)
    U = fit$U
    distr = dist(U)
    distc = dist(t(U))
    hr = hclust(distr)
    hc = hclust(distc)
    r = cutreeHybrid(hr, distM = as.matrix(distr),verbose=0)$labels
    c = cutreeHybrid(hc, distM = as.matrix(distc),verbose=0)$labels
    mats = rc.to.mats(X,r,c)
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  } 
  else if(method=="cobra"){
    n=nrow(X)
    p=ncol(X)
    A = B
    A = A - sum(A)/length(A)
    A = A/norm(A,"f")
    phi = 0.5; k = 5;
    
    library("cvxbiclustr")
    w_row = wts$w_row
    w_col = wts$w_col
    E_row = wts$E_row
    E_col = wts$E_col
    #### Initialize path parameters and structures
    nGamma = 30
    ## Generate solution path
    if(is.null(gamma)){
      gammaSeq <- 10**seq(.9,2,length.out=nGamma)
      sol = cobra_validate(A,E_row,E_col,w_row,w_col,gammaSeq)
      verr = sol$validation_error
      plot(verr)
      ix = which.min(verr)
      U = sol$U[[ix]]
    }
    else{
      sol = cobra_validate(A,E_row,E_col,w_row,w_col,gamma,fraction=0.01)
      U = sol$U[[1]]
    }
    sol$groups_row
    distr = dist(U)
    distc = dist(t(U))
    hr = hclust(distr)
    hc = hclust(distc)
    r = cutreeHybrid(hr, distM = as.matrix(distr),verbose=0)$labels
    c = cutreeHybrid(hc, distM = as.matrix(distc),verbose=0)$labels
    mats = rc.to.mats(X,r,c)
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  } else if(method=="spectral"){
    a = bicluster(X,as.integer(k))
    mats = rc.to.mats(X,a[[1]],a[[2]])
    res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  } else{
    print("method not supported")
  }
  res@RowxNumber[is.na(res@RowxNumber)]=FALSE
  res@NumberxCol[is.na(res@NumberxCol)]=FALSE
  ns = colSums(res@RowxNumber)
  res@RowxNumber = as.matrix(res@RowxNumber[,ns>1],ncol=sum(ns>1& ns < 200))
  res@NumberxCol = as.matrix(res@NumberxCol[ns>1,],ncol=sum(ns>1& ns < 200))
  if(nrow(res@NumberxCol)==1) res@NumberxCol=as.matrix(res@NumberxCol, ncol=1)
  res@Number=sum(ns>1 & ns < 200)
  res = BiclustResult(list(),res@RowxNumber,res@NumberxCol,res@Number,list())
  return(res)
}

bc.plot = function(X,res,labCol=NULL, labRow=NULL, rcols=NULL,colCol=NULL, colRow=NULL,mu=FALSE, gcols=NULL,...){
  B = res.to.B(X,res)
  Mu = res.to.Mu(X,res)
  if(is.null(rcols)) rcols = c("red","blue","green","orange","purple","cyan","brown","grey")
  if(is.null(gcols)){
    if(length(unique(labRow))>8) gcols = heat.colors(length(unique(labRow)))
    else gcols = rcols
  }
  c=order(Mu[1,])
  r=order(Mu[,1])
  coln = as.numeric(as.factor(Mu[1,]))
  rown = as.numeric(as.factor(Mu[,1]))
  if(mu){
    A=Mu
  } else{
    A=X
  }
  if(!is.null(labCol)){
    colm = as.numeric(as.factor(labCol))
    ColSideColors = rcols[colm][c]
    if(!is.null(labRow)){
      rowm = as.numeric(as.factor(labRow))
      RowSideColors = gcols[rowm][r]
      matshow(A[r,c], ColSideColors = ColSideColors, RowSideColors = RowSideColors,...)
    }
    else{
      matshow(A[r,c], ColSideColors = ColSideColors,...)
    }
  }
  else if(!is.null(labRow)){
    rowm = as.numeric(as.factor(labRow))
    RowSideColors = gcols[rowm][r]
    matshow(A[r,c], RowSideColors = RowSideColors,...)
  }
  else{
    matshow(A[r,c],...)
  }
}

matshow = function (A, key=FALSE, labRow = FALSE,labCol = FALSE,...){
  cols = rev(usecol(c(rev(  brewer.pal(n = 8, name = "Reds")), "white", brewer.pal(n = 8, name = "Blues")), n = 256))
  if(key){
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none",labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,colsep=c(0,ncol(A)),rowsep=c(0,nrow(X)),sepcolor="black",key=TRUE,...)
  }else{
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none", labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,  
              colsep=c(0,ncol(A)), rowsep=c(0,nrow(A)),sepcolor="black",keysize = .1,margins = c(2,2),key=FALSE,key.par = list(cex=.5),...)
  }
}

makekey = function (A, key=FALSE){
  z = seq(-max(abs(A)), max(abs(A)), length.out=100)
  breaks = seq(-max(abs(A)), max(A), length = 257)
  image(z = matrix(z, ncol = 1), col = cols, breaks = breaks, xaxt = "n", yaxt = "n")
  lv = pretty(breaks)
  minx = min(lv)
  ranx = diff(range(lv))
  xv = (lv - minx)/ranx
  xargs = list(at = xv, labels = lv, side=1)
  key.xlab = "Value"
  mtext(side = 1, key.xlab, line = par("mgp")[1],  padj = 0.5, cex = par("cex") * par("cex.lab"))
  title("")
  do.call(axis, xargs)
}


# Helper functions
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

matshow = function (A, key=FALSE, labRow = FALSE,labCol = FALSE,...){
  cols = rev(usecol(c(rev(  brewer.pal(n = 8, name = "Reds")), "white", brewer.pal(n = 8, name = "Blues")), n = 256))
  if(key){
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none",labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,colsep=c(0,ncol(A)),rowsep=c(0,nrow(X)),sepcolor="black",key=TRUE,...)
  } else{
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none", labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,  
              colsep=c(0,ncol(A)), rowsep=c(0,nrow(A)),sepcolor="black",keysize = .1,margins = c(2,2),key=FALSE,key.par = list(cex=.5),...)
  }
}

rc.to.mats = function(X, r, c){
  n = nrow(X)
  p = ncol(X)
  Mu = matrix(NA,n,p)
  for(i in 1:n){
    for(j in 1:p){
      Mu[i,j]=paste0(r[i],"_",c[j])
    }
  }
  Mu=matrix(as.numeric(as.factor(Mu)),n,p)
  K=min(c(max(Mu),50))
  rmat = matrix(FALSE,nrow=n,ncol=K)
  cmat = matrix(FALSE,nrow=K,ncol=p)
  for(k in 1:K){
    for(i in 1:nrow(Mu)){
      rmat[i,k] = as.logical(max(Mu[i,]==k))
    }
    for(j in 1:ncol(Mu)){
      cmat[k,j] = as.logical(max(Mu[,j]==k))
    }
  }
  return(list(rmat=rmat,cmat=cmat))
}

res.to.B = function(X, res){
  n = nrow(X)
  p = ncol(X)
  K = res@Number
  if(K==0) return(matrix(0, n, p))
  if(K==1) cmat=t(res@NumberxCol)
  else cmat = res@NumberxCol
  B = matrix(0,n,p)
  for(i in 1:n){
    for(j in 1:p){
      for(k in 1:K){
        if(res@RowxNumber[i,k] & cmat[k,j]) B[i,j] = k
      }
    }
  }
  return(B)
}

res.to.Mu = function(X, res){
  n = nrow(X)
  p = ncol(X)
  K = res@Number
  if(K==0) return(matrix(0, n, p))
  if(K==1) cmat=t(res@NumberxCol)
  B = res.to.B(X, res)
  Mu = matrix(0,n,p)
  for(k in 1:K){
    Mu[B==k]=mean(X[B==k])
  }
  
  return(Mu)
}

bc.plot = function(X,res,labCol=NULL, labRow=NULL, rcols=NULL,colCol=NULL, colRow=NULL,mu=FALSE, gcols=NULL,...){
  B = res.to.B(X,res)
  Mu = res.to.Mu(X,res)
  if(is.null(rcols)) rcols = c("red","blue","green","orange","purple","cyan","brown","grey")
  if(is.null(gcols)){
    if(length(unique(labRow))>8) gcols = heat.colors(length(unique(labRow)))
    else gcols = rcols
  }
  c=order(Mu[1,])
  r=order(Mu[,1])
  coln = as.numeric(as.factor(Mu[1,]))
  rown = as.numeric(as.factor(Mu[,1]))
  if(mu){
    A=Mu
  } else{
    A=X
  }
  if(!is.null(labCol)){
    colm = as.numeric(as.factor(labCol))
    ColSideColors = rcols[colm][c]
    if(!is.null(labRow)){
      rowm = as.numeric(as.factor(labRow))
      RowSideColors = gcols[rowm][r]
      matshow(A[r,c], ColSideColors = ColSideColors, RowSideColors = RowSideColors,...)
    }
    else{
      matshow(A[r,c], ColSideColors = ColSideColors,...)
    }
  }
  else if(!is.null(labRow)){
    rowm = as.numeric(as.factor(labRow))
    RowSideColors = gcols[rowm][r]
    matshow(A[r,c], RowSideColors = RowSideColors,...)
  }
  else{
    matshow(A[r,c],...)
  }
}

# Fit ssbc model
fit.ssbc = function(X,gamma1=.1,gamma2=.1,rho=.1,y=NULL,phi=.5,alpha=.5,k.row=10,k.col=10,eta=1,maxiter=10,tol.abs=1e-4,tol.rel=.01){
  
  # Parameter options
  supervised = !is.null(y)
  
  # objective function
  omegarow = function(U,wts){
    D = as.matrix(dist(U))
    E_row = wts$E_row
    w_row = wts$w_row
    s=0
    for(i in 1:nrow(E_row)){
      idx = which(as.logical(abs(E_row[i,])))
      s = s + w_row[i] * (D[idx[1],idx[2]])
    }
    return(s)
  }
  
  omegacol = function(U,wts){
    D = as.matrix(dist(t(U)))
    E_col = wts$E_col
    w_col = wts$w_col
    s=0
    for(i in 1:nrow(E_col)){
      idx = which(as.logical(abs(E_col[i,])))
      s = s + w_col[i] * (D[idx[1],idx[2]])
    }
    return(s)
  }
  
  obj = function(U,wts,beta){
    .5 * nrmF(X-U)^2 + gamma1 * (omegarow(U,wts) + omegacol(U,wts)) + gamma2 * sum(beta * apply(U,2,nrm2))
  }  
  
  U.update = function(U,X,V,Z,gamma2,beta,rho=1,niter=100){
    Ucol.update = function(u,x,v,z,gamma2,betai,rho=1,niter=100,tol=0.001){
      lambda = betai * gamma2
      f1 = function(u) .5*nrm2(x-u)^2 + (rho/2) * nrm2(u - v + z)^2
      f2 = function(u)  lambda * nrm2(u)
      df1 = function(u) (x-u) + rho*(u-v+z)
      prob = problem(f1=f1,f2=f2,df1=df1,proxf2=l2prox,lambda = lambda)
      sol = solve(prob,method="pgd",par=u,control=list(stepmethod="backtrack",maxiter=100,tol=1e-4))
      return(sol$x)
    }
    n = nrow(U)
    p = ncol(U)
    newU = matrix(NA,n,p)
    for(i in 1:p){
      newU[,i] = Ucol.update(U[,i],X[,i],V[,i],Z[,i],gamma2=gamma2,betai = beta[i],rho=rho,niter=niter)
    }
    return(newU)
  }
  
  V.update = function(U,Z,gamma1,rho,w_row,w_col,E_row,E_col){
    res = quiet(cobra(U+Z, E_row, E_col, w_row, w_col,gamma = gamma1/rho))
    return(res$U[[1]])
  }
  
  Z.update = function(U,V,Z){
    return(Z+U-V)
  }
  
  # Initialize
  n=nrow(X)
  p=ncol(X)
  phi=.05
  wts <- gkn_weights(as.matrix(X),phi=phi,k_row=5,k_col=5)
  w_row <- wts$w_row
  w_col <- wts$w_col
  E_row <- wts$E_row
  E_col <- wts$E_col
  U=X
  V=X
  Z=X
  
  # Get supervised weights
  if(supervised){
    distmat = as.matrix(dist(y))
    distmat[distmat==0]=1e-16
    a=rep(NA,length(w_row))
    for(i in 1:nrow(E_row)){
      idx = which(as.logical(abs(E_row[i,])))
      w_row[i] = alpha * w_row[i]  + (1-alpha) * (1 / (distmat[idx[1],idx[2]])^eta)
    }
    w_row = w_row * (1 / (sqrt(n)*sum(w_row)))
  }
  
  # Get adaptive weights
  res = quiet(cobra(X, E_row, E_col, w_row, w_col,gamma = 20))
  beta = 1/(apply(res$U[[1]],2,nrm2)^2)
  beta=rep(1,ncol(U))
  
  # History
  objx = obj(U,wts,beta)
  objs = objx
  primal.rs = NA
  dual.rs = NA
  
  # Just softthreshold
  if(gamma1==0){
    print("gamma1 = 0: just soft thresholding...")
    for(i in 1:ncol(U)){
      x = X[,i]
      u = U[,i]
      f1 = function(u) .5*nrm2(x-u)^2 
      f2 = function(u)  nrm2(u)
      df1 = function(u) x-u
      prob = problem(f1=f1,f2=f2,df1=df1,proxf2=l2prox,lambda = gamma2)
      sol = solve(prob,method="pgd",par=u,control=list(stepmethod="backtrack",maxiter=100,tol=1e-4))
      U[,i] = sol$x
    }
  }
  
  # Just bicluster
  else if(gamma2==0){
    print("gamma2 = 0: just convex biclustering...")
    res = quiet(cobra(X, E_row, E_col, w_row, w_col,gamma = gamma1))
    U = res$U[[1]]
  }
  
  # Sparse biclustering
  else{
    print("Sparse biclustering...")
    for(iter in 1:maxiter){
      oldV = V
      V = V.update(U,Z,gamma1,rho,w_row,w_col,E_row,E_col)
      U = U.update(U,X,V,Z,gamma2,beta=beta)
      Z = Z.update(U,V,Z)
      
      # Stopping criterion
      # https://arxiv.org/pdf/1704.06209.pdf
      n1 = nrow(U)*ncol(X)
      p1 = n1
      primal.e = sqrt(p1)*tol.abs + tol.rel*max(c(nrmF(U), nrmF(V)))
      dual.e = sqrt(n1)*tol.abs + tol.rel*nrmF(Z)
      
      objx = obj(U,wts,beta)
      primal.r = U - V
      dual.r = rho * (V - oldV)
      
      objs = c(objs, objx)
      primal.rs = c(primal.rs, nrm2(primal.r))
      dual.rs = c(dual.rs, nrm2(dual.r))
      
      if(nrm2(primal.r) < primal.e & nrm2(dual.r) < dual.e) break
      if(nrm2(primal.r) > 10*nrm2(dual.r)) rho = 2*rho
      if(nrm2(dual.r) > 10*nrm2(primal.r)) rho = rho/2
      
    }
  } 
  distr = dist(U)
  distc = dist(t(U))
  hr = hclust(distr)
  hc = hclust(distc)
  r = cutreeHybrid(hr, distM = as.matrix(distr),verbose=0)$labels
  c = cutreeHybrid(hc, distM = as.matrix(distc),verbose=0)$labels
  mats = rc.to.mats(X,r,c)
  res = BiclustResult(list(),mats$rmat,mats$cmat,ncol(mats$rmat),list())
  return(list(res=res,U=U,objs=objs,primal.rs=primal.rs,dual.rs=dual.rs))
}

spbc = function(X,k1,k2,lambda=0){
  sol = sparseBC(X,k1,k2,lambda)
  r = sol$Cs
  c = sol$Ds
  n = nrow(X)
  p = ncol(X)
  df1 = data.frame(subj=1:n, cs = r)
  df2 = data.frame(feat=1:p, cf = c)
  df = expand.grid(subj=df1$subj, feat=df2$feat)
  df = merge(df,df1, by="subj")
  df = merge(df,df2, by="feat")
  df = df[order(df$subj, df$feat),]
  df$cb = as.numeric(as.factor(paste(df$cs,df$cf)))
  df = df[,c("subj","feat","cs","cf","cb")]
  return(list(df, sol$mus))
}

spbc.kr = function (x, k, r, lambda, percent = 0.1, trace = FALSE) 
{
  miss <- sample(1:(nrow(x) * ncol(x)), nrow(x) * ncol(x), 
                 replace = FALSE)
  numberoftimes <- 1/percent
  allresults <- array(NA, dim = c(numberoftimes, length(k), 
                                  length(r)))
  Cs.init <- matrix(NA, nrow = nrow(x), ncol = length(k))
  for (i in 1:length(k)) {
    Cs.init[, i] <- kmeans(x, k[i], nstart = 20)$cluster
  }
  Ds.init <- matrix(NA, nrow = ncol(x), ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] <- kmeans(t(x), r[j], nstart = 20)$cluster
  }
  for (i in 1:numberoftimes) {
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss <- x
    missing <- miss[1:(nrow(x) * ncol(x))/numberoftimes]
    xmiss[missing] <- NA
    xmiss[missing] <- mean(xmiss, na.rm = TRUE)
    for (a in 1:length(k)) {
      for (b in 1:length(r)) {
        res <- sparseBC(xmiss, k[a], r[b], lambda = lambda, 
                        Cs.init = Cs.init[, a], Ds.init = Ds.init[, 
                                                                  b])$mus
        allresults[i, a, b] <- sum((x[missing] - res[missing])^2)
      }
    }
    miss <- miss[-1:-(dim(x)[1] * dim(x)[2]/numberoftimes)]
  }
  results.se <- apply(allresults, c(2, 3), sd)/sqrt(numberoftimes)
  results.mean <- apply(allresults, c(2, 3), mean)
  IndicatorMatrix <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                               1)] <= results.mean[2:length(k), 2:length(r)] + results.se[2:length(k), 
                                                                                                                          2:length(r)])
  if (max(IndicatorMatrix) == 0) 
    return(list(bestK = max(k), bestR = max(r)))
  RowIndexPlusColIndex <- outer(k[-length(k)], r[-length(r)], 
                                "*")
  smallestIndicatorTrue <- min(RowIndexPlusColIndex[IndicatorMatrix == 
                                                      TRUE])
  out <- which(IndicatorMatrix == TRUE & RowIndexPlusColIndex == 
                 smallestIndicatorTrue, arr.ind = TRUE)
  out <- matrix(c(k[out[, 1]], r[out[, 2]]), nrow(out), ncol(out))
  temprow <- NULL
  tempcol <- NULL
  for (i in 1:length(k)) {
    temprow <- c(temprow, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempcol <- c(tempcol, paste("R = ", r[i], sep = ""))
  }
  rownames(results.se) <- temprow
  colnames(results.se) <- tempcol
  rownames(results.mean) <- temprow
  colnames(results.mean) <- tempcol
  return(list(estimated_kr = out, results.se = results.se, 
              results.mean = results.mean))
}

spbc.cv = function(X,k1,k2,lambda=0){
  kr = spbc.kr(X,k1,k2,0,0.2)$estimated_kr
  sol = sparseBC(X,kr[1],kr[2],lambda)
  r = sol$Cs
  c = sol$Ds
  n = nrow(X)
  p = ncol(X)
  df1 = data.frame(subj=1:n, cs = r)
  df2 = data.frame(feat=1:p, cf = c)
  df = expand.grid(subj=df1$subj, feat=df2$feat)
  df = merge(df,df1, by="subj")
  df = merge(df,df2, by="feat")
  df = df[order(df$subj, df$feat),]
  df$cb = as.numeric(as.factor(paste(df$cs,df$cf)))
  df = df[,c("subj","feat","cs","cf","cb")]
  return(list(df, sol$mus))
}
