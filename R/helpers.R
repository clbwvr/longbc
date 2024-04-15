#' @export
matshow = function (A, key=FALSE, labRow = FALSE,labCol = FALSE,...){
  cols = rev(usecol(c(rev(brewer.pal(n = 7, name = "Reds")), "white", brewer.pal(n = 7, name = "Blues")), n = 256))
  #cols = usecol(c("#000089", "#0000ff","#8989ff","white", "#ff8989", "#ff0000","#890000"), n = 256)
  #cols = usecol(c("#0000ff","white",  "#ff0000"), n = 256)
  #cols = rbcols
  if(length(unique(A)) == 2) cols = rep("#FB7B5B",length(cols))
  if(key){
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none",labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,key=TRUE,...)
  }else{
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none", labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,
              sepcolor="black",keysize = .1,key=FALSE,key.par = list(cex=.5),...)
  }
}

f.plot=function(x,y,main="",...){
  par(tck=.03)
  plot(x, y,axes=FALSE,las=1,...)
  mtext(main, side=3, adj=0, line=1, cex = 2/3) 
  box(bty="L")
  par(family="Helvetica")
  y = c(y,2.2)
  py = pretty(y)
  px = pretty(x)
  axis(1,line=0,at=px,labels = px)
  axis(2,line=0,las=1,at=py,labels = py)
  dy = diff(py)[1]
  dx = diff(px)[1]
  rug(x = py + dy/2, ticksize = .015, side = 2,line=0)
  rug(x = px + dx/2, ticksize = .015, side = 1,line=0)
}

f.matplot=function(x,Y,main,ylim,...){
  par(tck=.03)
  par(mar=c(4,2.5,2.5,0.5),family="Helvetica")
  matplot(Y,las=1,axes=FALSE,ylim=ylim, ...)
  title(main, adj=0, line=1,font.main = 1)
  box(bty="L")
  par(family="Helvetica")
  #py = round(seq(-.6, .5,.2),1)
  py = pretty(seq(ylim[1], ylim[2],,100))
  px = 0:5
  axis(1,line=0,at=px+1,labels =0:5)
  axis(2,line=0,las=1,at=py,labels=py)
  dy = diff(py)[1]/2
  dx = diff(px)[1]/2
  rug(x = py + dy, ticksize = .015, side = 2,line=0)
  rug(x = px + dx, ticksize = .015, side = 1,line=0)
}

f.matplot1=function(x,Y,main,ylim,...){
  par(tck=.03)
  par(mar=c(4,4.5,2.5,2),family="Helvetica")
  matplot(Y,las=1,axes=FALSE,ylim=ylim, ...)
  title(main, adj=0, line=1,font.main = 1)
  box(bty="L")
  par(family="Helvetica")
  #py = round(seq(-.6, .5,.2),1)
  py = pretty(seq(ylim[1], ylim[2],,100))
  px = c(0,1,2)
  axis(1,line=0,at=c(1,2,3),srt=20,labels=c("Post-injury","Asymptomatic","Unrestricted"))
  axis(2,line=0,las=1,at=py,labels=py)
  dy = diff(py)[1]/2
  dx = diff(px)[1]/2
  rug(x = py + dy, ticksize = .015, side = 2,line=0)
  rug(x = px + dx, ticksize = .015, side = 1,line=0)
}

makekey = function (A){
  rrcols =  rev(usecol(c(rev(brewer.pal(n = 7, name = "Reds")), "white", brewer.pal(n = 7, name = "Blues")), n = 100))
  z = seq(-max(abs(A)), max(abs(A)), length.out=100)
  breaks = seq(-max(abs(A)), max(abs(A)), length = 101)
  image(z = matrix(z, nrow = 1), col = rrcols, breaks = breaks, xaxt = "n", yaxt = "n")
  lv = pretty(breaks)
  minx = min(lv)
  ranx = diff(range(lv))
  xv = (lv - minx)/ranx
  xargs = list(at = xv, labels = lv, side=4, las=1)
  key.xlab = ""
  title("")
  do.call(axis, xargs)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

logseq = function(from, to, length.out) exp(seq(log(from),log(to),length.out=length.out))

nrm1 = function(x) sum(abs(x))

nrm2 = function(x) as.numeric(sqrt(crossprod(x,x)))

l1prox = function(x,t){
  soft = function(u) sign(u)*max(abs(u)-t,0)
  return(sapply(x,soft))
}

l2prox = function(x,lam){
  nrm = nrm2(x)
  t = 1-lam/nrm
  s = max(0, t)
  return(s*x)
}

l12prox = function(x,lam1,lam2){
  l1prox(l2prox(x,lam2),lam1)
}

upper.to.sym = function(x){
  x[lower.tri(x)] = t(x)[lower.tri(x)]
  x
}
lower.to.sym = function(x){
  x[upper.tri(x)] = t(x)[upper.tri(x)]
  x
}

hush=function(code){
  sink("/dev/null")
  tmp = code
  sink()
  return(tmp)
}