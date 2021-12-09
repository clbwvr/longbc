context("longbc") 
testthat::test_that("longbc works and achieves ARI = 1 on simple example", {
  set.seed(123)
  parm = parmf(ns = 20, ps = 6, snrs = 2, nbs = 4, grids = TRUE,exchs = FALSE)
  datres = datf(parm[[1]])
  dat = datres[[2]]
  plotdat(dat)
  XtXrhoDeltai = datres[[3]]
  gtrue = dat[,c("subj","feat","cs","cf","cb")]
  gtrue = gtrue[!duplicated(gtrue),]
  lambdas = logseq(.01,1000,length.out=20)
  fit = longbc(dat = dat, Vstruct="id",lambdas=lambdas, tausd=1/4, rho=1, tol.abs=1e-2, tol.rel=1e-3,niter = 50, XtXrhoDeltai = XtXrhoDeltai,trace=TRUE,verbose=FALSE)
  best = ari.longbc(fit, gtrue)[[19]]
  testthat::expect_equal(best,1)
})