setwd("/users/cweaver3/projects/research/Longitudinal Biclustering/longbc/")
library(longbc)
source("sim_fns.R")
source("R/helpers.R")
logseq = function(from, to, length.out) exp(seq(log(from),log(to),length.out=length.out))
set.seed(123)
parm = parmf(ns = c(50,100), ps = c(10,30), snrs = c(1,2), nbs = c(9,16), ms=4, grids = c(TRUE, FALSE),exchs = c(TRUE, FALSE))
dats = mclapply(parm, datf,mc.cores = 7)

dat = datres[[1]]
XtXrhoDeltai = datres[[2]]
plotdat(dat)
gtrue = dat[,c("subj","feat","cs","cf","cb")]
gtrue = gtrue[!duplicated(gtrue),]
fit = longbc:::longbc(dat = dat, Vstruct="id",lambdas = seq(.1,5000,length.out=10), tausd=1/4, rho=1, tol.abs=1e-2, tol.rel=1e-3,niter = 50, XtXrhoDeltai = XtXrhoDeltai,trace=TRUE,nnk=5,wtype="l2",sparsew=FALSE,loud=FALSE)
longbc:::plot.longbc(fit,1)
longbc:::plot.longbc(fit,2)
longbc:::plot.longbc(fit,3)
longbc:::plot.longbc(fit,4)
longbc:::plot.longbc(fit,5,3)
longbc:::plot.longbc(fit,6,3)
which.max(longbc:::ari.longbc(fit, gtrue))
max(longbc:::ari.longbc(fit, gtrue))

dat1 =readRDS("../data/allqs_ts3_s.rds")
dat = dat1[,c("id","feature","tp","value")]
colnames(dat) = c("subj","feat","t","y")
dat$subj = as.numeric(factor(dat$subj))
dat$feat = as.numeric(factor(dat$feat))
dat$t = dat$t - 1
XtXrhoDeltai = readRDS("dtiXtXrhoDeltai.rds")
Rprof("a.out")
fit = longbc(dat = dat, Vstruct="id",lambdas = seq(0,400,length.out=8), XtXrhoDeltai = XtXrhoDeltai, tausd=1/2, rho=1, tol.abs=1e-1, tol.rel=1e-3,niter = 50, trace=TRUE,nnk=10,wtype="l2",sparsew=TRUE,loud=TRUE)
Rprof()
summaryRprof("a.out")
fit$dat$cb=1
plot(fit, 3)

