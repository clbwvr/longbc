setwd("/Users/cweaver3/Projects/Research/Longitudinal Biclustering/longbc")
source("bicluster_fns.r")
source("R/helpers.r")
library(longbc)
library(parallel)
library(lme4)
source("sim_fns.r")


parammat = readRDS("parammatl.rds")
parammat = parammat[parammat$exch==FALSE,]
paramlist = lapply(1:nrow(parammat), function(u) parammat[u,])
# saveRDS(paramlist, "paramlistl.rds")

# datf(paramlist[[1]])

### Simulate data
# simulated_data[ repetition ][ setting ][ param | data | XtXrhoDeltai | V | g]
# simulated_data = list()
# for(i in 1:20){
#   print(i)
#   simulated_data[[i]] = mclapply(paramlist, datf, mc.cores = detectCores()-1)
# }
# datf(paramlist[[30]])
# saveRDS(simulated_data, "simulated_data.rds")
simulated_data = readRDS("simulated_data.rds")

# ### Naive results
# # naive_results[ repetition ][ setting ][ data | nb | ari]
naive_results = list()
for(i in 1:length(simulated_data)){
  datsi = simulated_data[[i]]
  print(i)
  naive_results[[i]] = mclapply(datsi, function(u){
    tryCatch(dofnaive(u), error = function(e) NA)
  }, mc.cores = detectCores()-1)
}
dofnaive(simulated_data[[1]][[1]])
# saveRDS(naive_results, "naive_results.rds")
naive_results = readRDS("naive_results.rds")
#
b1=list()
for(i in 1:length(naive_results)){
  print(i)
  A = lapply(naive_results[[i]], function(u){
    tryCatch(return(u$ari), error=function(e) return(c(NA,NA,NA)))
  })
  b1[[i]] = do.call(rbind, A)
  b1[[i]] = cbind(parammat, b1[[i]])
  b1[[i]]$sim = i
}
aris = do.call(rbind, b1)
aris = aris[,colnames(aris)!="l"]
aris = aris[,colnames(aris)!="sim"]
m = melt(aris, id.vars = colnames(aris)[1:7])
ggplot(subset(m, variable!="1")) + geom_boxplot(aes(x = variable, y = value)) + facet_wrap(.~snr)
aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, length)
agg11=aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, mean)
agg12=aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, sd)
m1 = agg11
m2 = agg12
m1 = m1[order(m1$grid, m1$n, m1$p, m1$nb, m1$snr),]
m2 = m2[order(m2$grid, m2$n, m2$p, m2$nb, m2$snr),]
rounding = function(a) format(round(a, 2), nsmall = 2)
m3 = m1[,!(colnames(m1) %in% c("l","m","grid","exch"))]
m3$n = ifelse(m3$n != c(0, lag(m3$n)[-1]), m3$n,"")
m3$p = ifelse(m3$p != c(0, lag(m3$p)[-1]), m3$p,"")
m3$snr = ifelse(m3$snr != c(0, lag(m3$snr)[-1]), m3$snr,"")
m3$nb = ifelse(m3$nb != c(0, lag(m3$nb)[-1]), m3$nb,"")
m3$`1` = paste0(rounding(m1$`1`), " (", rounding(m2$`1`), ")")
m3$`2` = paste0(rounding(m1$`2`), " (", rounding(m2$`2`), ")")
m4 = m3[17:32,]
m3 = m3[1:16,]
m5 = cbind(m3, m4[,5:6])
m5 = m5[,c(1,2,4,3,5,6,7,8)]
kable(m5, booktabs = T, align = "c",row.names = F, format="latex")

### Proposed results
# prop_results[ repetition ][ setting ][ data | nb | ari]
prop_results = list()
for(i in 1:20){
  datsi = simulated_data[[i]]
  print(i)
  prop_results[[i]] = mclapply(1:4, function(u) dofprop(datsi[[u]], parammat[u,"l"]))
}
df = do.call(rbind, prop_results)
df = matrix(as.numeric(df),ncol=4)
colMeans(df)
saveRDS(prop_results, "prop_results.rds")
prop_results = readRDS("prop_results.rds")
b2=list()
for(i in 1:8){
  print(i)
  A = lapply(prop_results[[i]], function(u){
    if(is.null(names(u))){
      return(c(NA,NA,NA))
    } else{
      return(u$ari)
    }
  })
  b2[[i]] = do.call(rbind, A)
  b2[[i]] = cbind(parammat, b2[[i]])
  b2[[i]]$sim = i
}
aris = do.call(rbind, b2)
aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, length)
agg21 = aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, function(u) quantile(u, .7))
agg22 = aggregate(data=aris[,colnames(aris)!="sim"], . ~ n + p + snr + nb + m + grid + exch, sd)
m1 =  merge(agg11, agg21, by=c("n","p","snr","nb","m","grid","exch"))[,-c(8,12)]
m2 =  merge(agg12, agg22, by=c("n","p","snr","nb","m","grid","exch"))[,-c(8,12)]
m1 = agg11
m2 = agg12
rounding = function(a) format(round(a, 2), nsmall = 2)
colnames(m1)[9:10] = c("OLS","Long. BC")
colnames(m2)[9:10] = c("OLS","Long. BC")
m3 = m1
m3 = m3[order(m3$grid, m3$n, m3$p, m3$nb, m3$snr),]
m3 = m3[,c("n","p","nb","snr","OLS","Long. BC")]
m3$n = ifelse(m3$n != c(0, lag(m3$n)[-1]), m3$n,"")
m3$p = ifelse(m3$p != c(0, lag(m3$p)[-1]), m3$p,"")
m3$snr = ifelse(m3$snr != c(0, lag(m3$snr)[-1]), m3$snr,"")
m3$nb = ifelse(m3$nb != c(0, lag(m3$nb)[-1]), m3$nb,"")
m3$OLS = paste0(rounding(m1$OLS), " (", rounding(m2$OLS), ")")
m3$`Long. BC` = paste0(rounding(m1$`Long. BC`), " (", rounding(m2$`Long. BC`), ")")
m4 = m3[17:32,]
m3 = m3[1:16,]
m3$OLS1 = m4$OLS
m3$LongBC1 = m4$`Long. BC`

kable(m3, booktabs = T, align = "c",row.names = F, format="latex")
