library(pals)
library(gplots)
library(unikn)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(mvtnorm)
library(sparseBC)
library(ggplot2)
library(viridis)
source("/Users/cweaver3/Projects/Research/Longitudinal Biclustering/longbc/R/helpers.R")

# read data
setwd("/users/cweaver3/projects/Research/Longitudinal Biclustering/data")
allqs_ts = readRDS("allquantiles_ts.rds")
groups =  readRDS("groups.rds")

# only good tracts
goodlocs = c("unc","cgc","cgh","fma","fmi","cst","ar","atr","ptr","str","ilf","ifo","slf") 
allqs_ts$loc = sapply(allqs_ts$tract, function(u) strsplit(as.character(u),"_",TRUE)[[1]][[1]])
allqs_ts = subset(allqs_ts, loc %in% goodlocs)

# only subjects with time coverage
unq = allqs_ts[,c("id","tp")]
unq = unq[!duplicated(unq),]
tab = table(unq$id, unq$tp)
goodids = rownames(tab)[tab[,1]==1 & (tab[,2]==1 | tab[,3]==1)]
allqs_ts = subset(allqs_ts, id %in% goodids)
allqs_ts = allqs_ts[allqs_ts$id != "UCLA-FB-1111",]

# only good tracts with full coverage
allqs_ts$feature = paste0(allqs_ts$measure, ".", allqs_ts$tract)
unq = allqs_ts[,c("id","feature")]
unq = unq[!duplicated(unq),]
tab = table(unq$id, unq$feature)
goodfeats = colnames(tab)[colSums(tab) == max(colSums(tab))]
allqs_ts = subset(allqs_ts, feature %in% goodfeats)
p = length(unique(allqs_ts$feature))

# Invert FA
allqs_ts$value[allqs_ts$measure=="FA"] = 1 - allqs_ts$value[allqs_ts$measure=="FA"]

# Center at t = 3 and scale to SD=1
allqs_ts3 = subset(allqs_ts, tp<4)
allqs_ts3$feature = paste0(allqs_ts3$measure, ".", allqs_ts3$tract)
agg = aggregate(data=subset(allqs_ts3,tp==3), value~feature, mean)
colnames(agg)[2] = "featmean"
allqs_ts3_s = merge(allqs_ts3, agg, by="feature")
allqs_ts3_s$value = allqs_ts3_s$value-allqs_ts3_s$featmean
agg = aggregate(data=allqs_ts3_s, value~feature, sd)
colnames(agg)[2] = "featsd"
allqs_ts3_s = merge(allqs_ts3_s, agg, by="feature")
allqs_ts3_s$value = allqs_ts3_s$value/allqs_ts3_s$featsd
saveRDS(allqs_ts3_s,"allqs_ts3_s.rds")

# Empirical covariance
D = dcast(data=allqs_ts3_s, id+tp ~ feature)
D = data.frame(D)[,-c(1,2)]
D = D[complete.cases(D),]
Vcov = cov(D)
Vcor = cor(D)
tract = sapply(colnames(Vcov), function(u) strsplit(as.character(u),".",TRUE)[[1]][2])
meas = sapply(colnames(Vcov), function(u) strsplit(as.character(u),".",TRUE)[[1]][1])
cols = meas
cols = ifelse(cols=="Da",bcols[1], ifelse(cols=="Dr", bcols[2], ifelse(cols=="FA",bcols[3],bcols[4])))
cols = cols[order(meas)]
Vcov = Vcov[order(meas, tract),order(meas, tract)]
Vcor = Vcor[order(meas, tract),order(meas, tract)]
matshow(Vcor,ColSideColors= cols, RowSideColors=cols,rowsep=c(0,nrow(Vcov)),colsep=c(0,ncol(Vcov)))

setwd("../figs")
textwidth.dissertation  = 6.25;
pointsize.dissertation = 8
pdf("cor.pdf",.4 * textwidth.dissertation, .4 * textwidth.dissertation, pointsize = pointsize.dissertation)
par(mar=c(0,0,0,0),family="Helvetica")
matshow(Vcor,ColSideColors= cols, RowSideColors=cols)
dev.off()
pdf("keycor1.pdf", .10 * textwidth.dissertation, .4 * textwidth.dissertation, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(Vcor)
dev.off()

# Estimate covariance
ids = unique(allqs_ts_s$id)
res = list()
allsub = c()
n = length(ids)
for(i in 1:length(ids)){
  idi = ids[i]
  sub = subset(allqs_ts_s, id==idi)
  fit = lm(data=sub, value ~ feature:tp)
  sub$res = residuals(fit)
  allsub=rbind(sub,allsub)
}
makeP = function(theta){
  theta0 = theta[1]
  theta1 = theta[2]
  theta2 = theta[3]
  P = matrix(theta1, length(me),length(me))
  thetadiff = theta1 + theta2
  P[me=="Da",] = P[,me=="Da"] = thetadiff
  diag(P) = 1
  Pcov = theta0 * P
  Pcor = P
  return(list(Pcov, Pcor))
}

f = function(theta){
  P = makeP(theta)[[1]]
  -sum(dmvnorm(D,sigma = P,log=TRUE))
}
o = optim(c(1,0,0),f,control = list(maxit=1000))
theta = o$par
theta0 = theta[1]
theta1 = theta[2]
theta2 = theta[3]
Pcov = makeP(theta)[[1]]
Pcor = makeP(theta)[[2]]

pdf("estcor.pdf",.4 * textwidth.dissertation, .4 * textwidth.dissertation, pointsize = pointsize.dissertation)
par(mar=c(0,0,0,0),family="Helvetica")
matshow(Pcor,ColSideColors= cols, RowSideColors=cols)
dev.off()
pdf("keycor2.pdf", .10 * textwidth.dissertation, .4 * textwidth.dissertation, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(Pcor)
dev.off()

V1 = mean(diag(V))
Vo = V; diag(Vo) = NA;Vo[me!="Da",]=NA;
matshow(Vo)
Vda = mean(Vo[upper.tri(Vo)],na.rm=TRUE)
Vo = V; diag(Vo) = NA;Vo[me=="Da",]=NA;Vo[,me=="Da"]=NA;
Vdb = mean(Vo, na.rm=TRUE)
Vest = matrix(Vdb,length(me),length(me))
Vest[me=="Da",] = Vda
Vest[,me=="Da"] = Vda
diag(Vest) = V1
Vp = Vest / V1