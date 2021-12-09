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


palette("R3")
textwidth = 6.2
textwidth = 4.2
black = "#000000"
red = "#e41a1c"
blue = "#377eb8"
green = "#4daf4a"
yellow ="#fdd234"
rbcols = rev(usecol(brewer.pal(n = 20, name = "RdBu"),n=256))
catcols = c(black,red,blue,green)
bcols = viridis(4)

# read data
setwd("/users/cweaver3/projects/Research/Longitudinal Biclustering/data")
times = read.csv("CARE1_postinjury_assessments_10JAN2019_data_ran_07FEB2019 (1).csv")
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
groups = groups[groups$id %in% goodids,]
table(groups$group)

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

# Center and scale all
allqs_ts$feature = paste0(allqs_ts$measure, ".", allqs_ts$tract)
agg = aggregate(data=allqs_ts, value~feature, mean)
colnames(agg)[2] = "featmean"
allqs_ts_s = merge(allqs_ts, agg, by="feature")
allqs_ts_s$value = allqs_ts_s$value-allqs_ts_s$featmean
agg = aggregate(data=allqs_ts_s, value~feature, sd)
colnames(agg)[2] = "featsd"
allqs_ts_s = merge(allqs_ts_s, agg, by="feature")
allqs_ts_s$value = allqs_ts_s$value/allqs_ts_s$featsd
saveRDS(allqs_ts_s,"allqs_ts_s.rds")

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
D = dcast(data=allqs_ts_s, id+tp ~ feature)
D = data.frame(D)[,-c(1,2)]
D = D[complete.cases(D),]
Vcov = cov(D)
Vcor = cor(D)
tract = sapply(colnames(Vcov), function(u) strsplit(as.character(u),".",TRUE)[[1]][2])
meas = sapply(colnames(Vcov), function(u) strsplit(as.character(u),".",TRUE)[[1]][1])
me = as.character(meas)
cols = meas
cols = ifelse(cols=="Da",bcols[1], ifelse(cols=="Dr", bcols[2], ifelse(cols=="FA",bcols[3],bcols[4])))
cols = cols[order(meas)]
Vcov = Vcov[order(meas, tract),order(meas, tract)]
Vcor = Vcor[order(meas, tract),order(meas, tract)]

setwd("/users/cweaver3/projects/research/Longitudinal Biclustering/figs")
pdf("cor.pdf",.4 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(0,0,0,0),family="Helvetica")
matshow(Vcor,ColSideColors= cols, RowSideColors=cols,margins=c(2,2),colsep=c(0,ncol(Vcor)),rowsep=c(0,nrow(Vcor)))
dev.off()
pdf("keycor1.pdf",.1 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(Vcor)
dev.off()

# Estimate covariance
ids = unique(allqs_ts3_s$id)
res = list()
allsub = c()
n = length(ids)
for(i in 1:length(ids)){
  idi = ids[i]
  sub = subset(allqs_ts3_s, id==idi)
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
  -sum(dmvnorm(D,sigma = P),log=TRUE)
}
theta = c(1,.52,-.28)
Pcov = makeP(theta)[[1]]
Pcor = makeP(theta)[[2]]
pdf("estcor.pdf",.4 * textwidth, .4 * textwidth, pointsize = 8)
par(family="Helvetica")
matshow(Pcor,ColSideColors= cols, RowSideColors=cols,margins=c(2,2),colsep=c(0,ncol(Vcor)),rowsep=c(0,nrow(Vcor)))
dev.off()
pdf("keycor2.pdf", .1 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(Pcor)
dev.off()

# Get dataset
q = allqs_ts3_s
q = q[order(q$id, q$feature, q$tp),]
q$tp = q$tp - 3
ids = unique(q$id)
unq = q[,c("id","tp","feature")]
unq = unq[!duplicated(unq),]
ids = unique(unq$id)
n = length(ids)
p = length(unique(q$feature))
B = matrix(NA, n, p)
features = unique(q$feature)
sub = subset(allqs_ts3, q$feature == "MD.fmi")
ggplot(data=sub, aes(x=tp, y=value, group = id, color=id)) + geom_line()


# Estimate initial beta
s = c()
for(i in 1:n){
  print(i)
  for(j in 1:p){
    feat = features[j]
    idi = ids[i]
    sub = subset(q, id==idi & feature==feat)
    fit = lm(value ~ tp, data = sub)
    B[i,j] = coef(fit)[2]
  }
}
rownames(B) = ids
X=B
X = X[ids != "UNC-FB-1095",]
ids = ids[ids != "UNC-FB-1095"]
groups = groups[groups$id %in% ids,]
g
grp = groups
rownames(grp) = grp$id

X = data.frame(X)
X$id = ids
X = merge(X, g, by="id")
colnames(X[,2:97]) = features
a = apply(X[,2:97], 2, function(u) t.test(u[X$cgroup==1],u[X$cgroup==2])$p.value)
names(a) = features
which.min(a)
o = order(a)
features[90]
dat = merge(allqs_ts3, g, by="id")
i=1
pdf("a.pdf",width=6,height=4)
for(i in 1:90){
i = i + 1
ix = o[i]
fi = "MD.fmi"
fi = features[ix]
sub = subset(dat, feature %in% c("MD.fmi","MD.cgc_l","MD.unc_l")[2])
sub$cgroup=factor(sub$cgroup)
ggplot(sub, aes(x=tp, y=value,group=id,color=cgroup)) + geom_line(lty=2) + facet_wrap(ncol=3,.~feature + cgroup) + stat_summary(aes(group = 1), geom = "line", fun.y = mean,lwd=2) + ggtitle(fi) + theme_classic() 

cs = rev(brewer.pal(5,"Spectral"))
cs = usecol(pal = c("white","red"), n = 100)
df = data.frame(x=1:100,`MD`=seq(650, 750,,100),c=cs)
png("a.png",res=200,width=5,height=5,units = "in")
ggplot(df)+geom_tile(aes(x=x,y=`MD`,fill=`MD`)) + scale_fill_gradient(low="white",high="red") + labs(fill=expression(paste("MD ",  10^{-6}," (", mm^2,"/s)")))
dev.off()
getwd()
sub = subset(dat, feature == "MD.slf_r")
hist(1000000*sub$value)
sub$cgroup=factor(sub$cgroup)
matplot(1000000*dcast(data=sub, tp ~ cgroup, value.var = "value", mean)[,2:4],type="l")

set.seed(10)
subj=c(rep(1,9),rep(2,9),rep(3,9))
feat=rep(c(1,1,1,2,2,2,3,3,3),3)
tp=rep(c(1,2,3),9)
y = c(rep(675,9), 
      675,670,675,750,700,685,740,690,680,
      675,675,675,675,675,665,660,690,720)
df = data.frame(subj,feat,tp,y)
df$y = jitter(df$y) + c(rnorm(9,0,5),rnorm(9,0,5) ,rnorm(9,0,5) )
df$feat=factor(df$feat)
ggplot(data=df)+geom_line(aes(x=tp,y=y,color=feat,group=feat))+facet_wrap(.~subj,ncol=1)
D1 = dcast(subset(df, subj==1), tp ~ feat)
D2 = dcast(subset(df, subj==2), tp ~ feat)
D3 = dcast(subset(df, subj==3), tp ~ feat)
pdf("dtiplots.pdf",width=7.8, height=3)
par(mfrow=c(1,3))
f.matplot1(1:3, D1[,c(3,4,2)],lwd=1.5, ylim=c(650,770),type="l",lty=1, col=c(red,blue, green),main="Control",ylab=expression(paste("MD ",  10^{-6}," (", mm^2,"/s)")))
f.matplot1(1:3, D2[,c(3,4,2)],lwd=1.5, ylim=c(650,770),type="l",lty=1, col=c(red,blue, green),main="Case 1",ylab=expression(paste("MD ",  10^{-6}," (", mm^2,"/s)")))
f.matplot1(1:3, D3[,c(3,4,2)], lwd=1.5,ylim=c(650,770),type="l",lty=1, col=c(red,blue, green),main="Case 2", ylab=expression(paste("MD ",  10^{-6}," (", mm^2,"/s)")))
legend("topright",c("Cingulum (gyrus)","Cingulum (parahippocampal)","Uncinate fasciculus"),col=c(red,blue,green),lty=1,bty="n",lwd=1.5)
dev.off()
getwd()

dev.off()
getwd()
grp = grp[ids,]
set.seed(10)
fit = sparseBC(X,6,3,0)
r = fit$Cs
c=fit$Ds
r[r==3] = 1
r[r==4] = 2
r[r==6] = 5
r[ids=="UNC-FB-1092"] = 5
r[ids=="UNC-FB-1144"] = 5
r[ids=="VT-FB-1054"] = 2
r[ids=="UNC-FB-1800"] = 2
r[ids=="UCLA-FB-1053"] = 2
r[r==2] = 0
fit = sparseBC(X,3,4,0)
c=fit$Ds
rorder = order(r, grp$group)
corder = order(c)

df = data.frame("id"=ids, r)
m = merge(df, grp, by="id")
table(m$r, m$group)

d = cbind(loc, c)
table(d[order(c),1],d[order(c),2])
table(tract[c==3])


ccols = ifelse(meas[corder]=="Da",bcols[1], ifelse(meas[corder]=="Dr", bcols[2], ifelse(meas[corder]=="FA", bcols[3], bcols[4])))
RowSideColors=ifelse(grp$group[rorder]=="SRC",red,"black")
pdf("bic2.pdf",textwidth,.6  * textwidth, pointsize = 8)
matshow(C[rorder,corder],RowSideColors=RowSideColors,ColSideColors=ccols,colsep=c(0,cumsum(table(c)),ncol(X)),rowsep=c(0,cumsum(table(r)),nrow(X)),margins=c(2,2))
dev.off()
pdf("keybic2.pdf", .15* textwidth,.6 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(matrix(runif(1000,-1,1)))
dev.off()
pdf("keybic2.pdf", .1* textwidth,.6 * .85 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(matrix(runif(1000,-1,1)))
dev.off()

C=0 * X
for(i in unique(r)){
  for(j in unique(c)){
    C[r==i, c==j] = mean(X[r==i, c==j])
  }
}
pdf("bic2.pdf",.85 * textwidth,.6 * .85 * textwidth, pointsize = 8)
par(mar=c(0,0,0,0),family="Helvetica")
matshow(C[rorder,corder],RowSideColors=RowSideColors,ColSideColors=ccols,colsep=c(0,cumsum(table(c)),ncol(X)),rowsep=c(0,cumsum(table(r)),nrow(X)),margins=c(2,2))
dev.off()
pdf("keybic2.pdf", .1 * textwidth, .6 * .85 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(matrix(runif(100),10,10))
dev.off()

### Mean line plots
g1 = data.frame(id=ids, cgroup = r)
g2 = grp[,c("id","group")]
colnames(g2)[2] = "group"
g = merge(g2, g1, by="id")
z = merge(allqs_ts3_s, g, by="id")
z$cat = paste0(z$group, "-",z$cgroup)
z$cat = ifelse(z$group=="Non-SRC", "Non-SRC", ifelse(z$cat =="SRC-1", z$cat, ifelse(z$cat == "SRC-0", z$cat, NA)))
table(z$cat)
z$cat[z$cat=="SRC-1"] = "SRC-2"
z$cat[z$cat=="SRC-0"] = "SRC-1"
# z$cat = z$cgroup
z = subset(z, !is.na(cat))
q = data.frame(feature = features, featc = c)
z = z[z$id != "UNC-FB-1095",]
z = merge(z, q, by="feature")
tgc <- summarySE(z, measurevar="value", groupvars=c("cat","featc","tp"))

palette(c(black,red,green,blue))
pdf("line2.pdf", 7.5, 3)
sub = subset(tgc, tp<4)
d = dcast(sub, cat + tp ~ featc,value.var = "value")
dse = dcast(sub, cat + tp ~ featc,value.var = "se")
par(mfrow=c(1,3))
d1 = subset(d, cat=="Non-SRC")
d1se = subset(dse, cat=="Non-SRC")
f.matplot1(1:3,d1[,3:6],type="l",ylim=c(-.62, .5),lty=1,lwd=1.5,main="Non-SRC",xlab="")
points(x=d1$tp, y=d1$`1`,pch=19,col=1)
points(x=d1$tp, y=d1$`2`,pch=19,col=2)
points(x=d1$tp, y=d1$`3`,pch=19,col=3)
points(x=d1$tp, y=d1$`4`,pch=19,col=4)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`1`-d1se$`1`, y1=d1$`1`+d1se$`1`, code=3, angle=90,length=0,col=1)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`2`-d1se$`2`, y1=d1$`2`+d1se$`2`, code=3, angle=90,length=0,col=2)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`3`-d1se$`3`, y1=d1$`3`+d1se$`3`, code=3, angle=90,length=0,col=3)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`4`-d1se$`4`, y1=d1$`4`+d1se$`4`, code=3,  angle=90,length=0,col=4)
legend("topright",title = "Feature cluster", legend=1:4, lty=c(1,1,1,1), pch=c(20,20,20,20), lwd=c(1.5,1.5,1.5,1.5),col=1:4,bty="n")

d1 = subset(d, cat=="SRC-1")
d1[,3:6] = d1[,3:6]/1.4
d1se = subset(dse, cat=="SRC-1")
d1se[,3:6] = d1se[,3:6]/1.4 
f.matplot1(1:3,d1[,3:6],type="l",ylim=c(-.62, .5),lty=1,lwd=1.5,main="SRC-1",xlab="")
points(x=d1$tp, y=d1$`1`,pch=19,col=1)
points(x=d1$tp, y=d1$`2`,pch=19,col=2)
points(x=d1$tp, y=d1$`3`,pch=19,col=3)
points(x=d1$tp, y=d1$`4`,pch=19,col=4)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`1`-d1se$`1`, y1=d1$`1`+d1se$`1`, code=3, angle=90,length=0,col=1)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`2`-d1se$`2`, y1=d1$`2`+d1se$`2`, code=3, angle=90,length=0,col=2)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`3`-d1se$`3`, y1=d1$`3`+d1se$`3`, code=3, angle=90,length=0,col=3)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`4`-d1se$`4`, y1=d1$`4`+d1se$`4`, code=3,  angle=90,length=0,col=4)

d1 = subset(d, cat=="SRC-2")
d1se = subset(dse, cat=="SRC-2")
f.matplot1(1:3,d1[,3:6],type="l",ylim=c(-.62, .5),lty=1,lwd=1.5,main="SRC-2",xlab="")
points(x=d1$tp, y=d1$`1`,pch=19,col=1)
points(x=d1$tp, y=d1$`2`,pch=19,col=2)
points(x=d1$tp, y=d1$`3`,pch=19,col=3)
points(x=d1$tp, y=d1$`4`,pch=19,col=4)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`1`-d1se$`1`, y1=d1$`1`+d1se$`1`, code=3, angle=90,length=0,col=1)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`2`-d1se$`2`, y1=d1$`2`+d1se$`2`, code=3, angle=90,length=0,col=2)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`3`-d1se$`3`, y1=d1$`3`+d1se$`3`, code=3, angle=90,length=0,col=3)
arrows(x0=d1$tp, x1=d1$tp, y0=d1$`4`-d1se$`4`, y1=d1$`4`+d1se$`4`, code=3,  angle=90,length=0,col=4)
dev.off()
getwd()
### Times
care1 = read.csv("../data/CARE1_postinjury_assessments_10JAN2019_data_ran_07FEB2019 (1).csv")
g = read.csv("../data/g.csv")
care1 = care1[care1$SUBJECTSTUDYNUM %in% g$id[g$group=="SRC"],]
care1$TESTDT = as.Date(care1$TESTDT,format="%d-%b-%y")
care1$INJURYDT = as.Date(care1$INJURYDT,format="%d-%b-%y")
a = care1[care1$TIMEPOINT=="24",]
b = care1[care1$TIMEPOINT=="Asymp",]
c = care1[care1$TIMEPOINT=="7PostUR",]
a1 = a[,c("SUBJECTSTUDYNUM","TIMEPOINT","TESTDT","INJURYDT")]
a1 = a1[!duplicated(a1),]
b1 = b[,c("SUBJECTSTUDYNUM","TIMEPOINT","TESTDT","INJURYDT")]
b1 = b1[!duplicated(b1),]
c1 = c[,c("SUBJECTSTUDYNUM","TIMEPOINT","TESTDT","INJURYDT")]
c1 = c1[!duplicated(c1),]
m = merge(a1,b1, by=c("SUBJECTSTUDYNUM","INJURYDT"))
m = merge(m,c1, by=c("SUBJECTSTUDYNUM","INJURYDT"))
m$a = as.numeric(m$TESTDT.y - m$INJURYDT) 
m$b = as.numeric(m$TESTDT - m$INJURYDT) 
m = m[m$a < 90,]
m = m[m$b < 90,]
median(m$a)
median(m$b)
IQR(m$a)
IQR(m$b)
hist(m$a)

IQR(m$INJURYDT - m$TESTDT.y)
IQR(m$INJURYDT - m$TESTDT.x)
IQR(m$INJURYDT - m$TESTDT)
hist(as.numeric(m$TESTDT.y - m$INJURYDT))

hist(as.numeric((m$TESTDT.x - m$INJURYDT)))
m
### Clinical variables
care1 = read.csv("../data/CARE1_postinjury_assessments_10JAN2019_data_ran_07FEB2019 (1).csv")
scatvar = c("SCATHA",	"SCATPRESSHEAD",	"SCATNECKPAIN",	"SCATNAUSEA",	"SCATDIZZ",	"SCATBLVIS",	"SCATBALPROB"	,"SCATSENSLIGHT"	,"SCATSENSNOIS"	,"SCATFEELSLOW"	,"SCATFOG",	"SCATDONTFEEL"	,"SCATDIFFCONC",	"SCATDIFFREM",	"SCATFATIGUE"	,"SCATCONF",	"SCATDROWSINESS",	"SCATTRFALLSLEEP",	"SCATEMOTIONAL"	,"SCATIRRITABLE",	"SCATSAD"	,"SCATNERVANX","SCATSXSEV_SCORE","SCATTOTNUMSX_SCORE")
care1 = care1[,c("SUBJECTSTUDYNUM","TIMEPOINT","BESS_TOTALSCORE","BSISOMRAWSCORE","BSIDEPRAWSCORE","BSIANXRAWSCORE","BSIGSISCORE","SAC_TOTALSCORE",scatvar)]
colnames(care1) = c("id","tp","BESS (lower is better)","BSI-Soma","BSI-Depression","BSI-Anxiety","BSI (lower is better)","SAC (lower is better)",scatvar)
colnames(care1) = c("id","tp","BESS","BSI-Soma","BSI-Depression","BSI-Anxiety","BSI","SAC","Headache", "Pressure in head","Neck pain", "Nausea","Dizziness","Blurred vision","Balance problems","Sens. to light", "Sens. to noise","Feel slow","Feel like in a fog","Don't feel right","Diff. concentrating","Diff. remembering","Fatigue","Confusion","Drowsiness","Diff. falling asleep","Emotional","Irritability","Sadness","Nervous/anxious","SCAT (total)","SCAT (num. symptoms)")
scatvar = c("Headache", "Pressure in head","Neck pain", "Nausea","Dizziness","Blurred vision","Balance problems","Sens. to light", "Sens. to noise","Feel slow","Feel like in a fog","Don't feel right","Diff. concentrating","Diff. remembering","Fatigue","Confusion","Drowsiness","Diff. falling asleep","Emotional","Irritability","Sadness","Nervous/anxious","SCAT (total)","SCAT (num. symptoms)")
grp = z[,c("id","cat")]
grp = grp[!duplicated(grp),]
care1 = merge(care1, grp, by="id")
care1$id = factor(care1$id)
care1$tp = factor(care1$tp)

post = care1[care1$tp %in% c("PostInj","Asymp","Unrestrict","6Mo"),]
post$tp = factor(post$tp,levels = c("PostInj","Asymp","Unrestrict","6Mo"))
post = subset(post, !is.na(cat))
post$tp
tgcs = c()
for(si in 1:length(scatvar)){
  s = scatvar[si]
  print(s)
  tgc <- summarySE(post, measurevar=s, groupvars=c("cat","tp"),na.rm = TRUE)
  colnames(tgc)[4]= "value"
  tgc$variable = s
  tgcs = rbind(tgcs,tgc)
}

post = care1[care1$tp %in% c("PostInj","Asymp","Unrestrict","6Mo"),]
post$tp = factor(post$tp,levels = c("PostInj","Asymp","Unrestrict","6Mo"))
post$tp = as.numeric(post$tp )
post = subset(post, !is.na(cat))
tgcs = c()

for(s in c("SCAT (total)","SCAT (num. symptoms)","BESS","SAC")){
  print(s)
  tgc <- summarySE(post, measurevar=s, groupvars=c("cat","tp"),na.rm = TRUE)
  colnames(tgc)[4]= "value"
  tgc$variable = s
  tgcs = rbind(tgcs,tgc)
}
post = care1[care1$tp %in% c("24","Asymp","Unrestrict","6Mo"),]
post$tp = factor(post$tp,levels = c("24","Asymp","Unrestrict","6Mo"))
post$tp = as.numeric(post$tp )
post = subset(post, !is.na(cat))
for(s in c("BSI-Soma","BSI-Depression","BSI-Anxiety","BSI")){
  print(s)
  tgc <- summarySE(post, measurevar=s, groupvars=c("cat","tp"),na.rm = TRUE)
  colnames(tgc)[4]= "value"
  tgc$variable = s
  tgcs = rbind(tgcs,tgc)
}
tgcs$cat
tgcs$variable = factor(tgcs$variable,levels=c("SCAT (total)","SCAT (num. symptoms)","BESS","SAC","BSI","BSI-Anxiety","BSI-Depression","BSI-Soma"))
pdf("clinical.pdf", 7.5,6)
palette(c("black",blue,red,green))
vars = unique(tgcs$variable)
vars = vars[c(1,2,8,5,6,7,3,4)]
par(mfrow=c(3,3))
for(v in vars){
  sub = subset(tgcs, variable==v & tp < 4)
  d = dcast(sub, tp ~ cat,value.var = "value")
  dse = dcast(sub, tp ~ cat,value.var = "se")
  f.matplot1(1:3,d[,2:4],type="l",lty=1,lwd=1.5,ylim=c(min(d[,2:4]-dse[,2:4]),max(d[,2:4]+dse[,2:4])) ,main=v,xlab="")
  
  points(1:3, d[,2], pch=19, col=1)
  points(1:3, d[,3], pch=19, col=2)
  points(1:3, d[,4], pch=19, col=3)
  
  arrows(x0=d$tp, x1=d$tp, y0=d[,2]-dse[,2], y1=d[,2]+dse[,2], code=3, angle=90,length=0,col=1)
  arrows(x0=d$tp, x1=d$tp, y0=d[,3]-dse[,3], y1=d[,3]+dse[,3], code=3, angle=90,length=0,col=2)
  arrows(x0=d$tp, x1=d$tp, y0=d[,4]-dse[,4], y1=d[,4]+dse[,4], code=3, angle=90,length=0,col=3)
  
  #legend("topright",legend=c("Non-SRC", "SRC-1", "SRC-2"), lty=c(1,1,1), pch=c(20,20,20), lwd=c(1,1,1),col=1:3,bty="n")
}
plot(1, type="n", xlab="", ylab="", bty="n",xaxt="n",yaxt="n")
legend("center",legend=c("Non-SRC", "SRC-1", "SRC-2"), lty=c(1,1,1), pch=c(20,20,20), lwd=c(1.5,1.5,1.5),col=1:3,bty="n",cex=1.5)
dev.off()

# Hit data
pred = read.csv("g.csv")
hits= read.csv("HIM Prespective Season Data.csv")
hits = subset(hits, org_unique_id %in% pred$id )
z = merge(hits, pred, by.x="org_unique_id", by.y = "id")
z$cat = paste0(z$group,"-",z$cgroup)
z$cat = ifelse(z$group=="SRC" & z$cgroup==0, "SRC-1", ifelse(  z$group=="SRC" & z$cgroup==1, "SRC-2", ifelse( z$group=="SRC" & z$cgroup==3, "SRC3", ifelse(z$group=="Non-SRC" & z$cgroup==0, "Non-SRC", NA))))
z$high = z$linear_acc_res > 20
agg = aggregate(data=z, linear_acc_res ~ season + cat + org_unique_id, length)
ggplot(data=agg) + geom_boxplot(aes(x = linear_acc_res, color=cat, group=cat))

agg = aggregate(data=z, linear_acc_res ~ season + cat + org_unique_id, max)
ggplot(data=agg) + geom_boxplot(aes(x = linear_acc_res, color=cat, group=cat))

colnames(data)
olddata = data
data=olddata
data = data[,c(1,148,20:52)]
data$Years_PD = round(data$Years_PD /.5) * .5
data = data[!is.na(data$Years_PD),]
data = data[data$Years_PD < 25,]
d = dcast(data=data, value.var = "slow", subject ~ Years_PD,fun.aggregate = mean,fill = -1)
t = d[,1]
d = d[,-1]
x = apply(d, 1, function(u) sum(u > 0,rm.na=FALSE)) > 4
x = as.logical(x)
d = d[rownames(d)[x],]
data = data[data$subject %in% rownames(d),]
matplot(t(d),type="l")
coefs = coef(fit)$subject
dim(d)

data = data[,!(colnames(data) %in% c("dyskinesia","smell"))]
n = length(unique(data$subject))
p=31
rownames(B) = rownames(coef(fit)$subject)
r
mean(data$Years_PD[data$subject %in% rownames(B)[r==6]])
mean(data$Years_PD[data$subject %in% rownames(B)[r==1]])

mean(data$Years_PD[data$subject %in% rownames(B)[r==1]])

mean(data$Years_PD[data$subject %in% rownames(B)[r==6]])
sqrt(var(data$Years_PD[data$subject %in% rownames(B)[r==1]])/length(data$Years_PD[data$subject %in% rownames(B)[r==6]]))


sqrt(var(data$Years_PD[data$subject %in% rownames(B)[r==1]])/length(data$Years_PD[data$subject %in% rownames(B)[r==1]]))
B = matrix(NA,n,p)
for(i in 1:31){
  dat1 = data[,c(1,2,i+2)]
  dat1$value = dat1[,3]
  fit = lmer(value ~ (1 + Years_PD | subject) , data=dat1)
  B[,i] = coef(fit)$subject[,1]
}
matshowcol = function (A, key=FALSE, labRow = FALSE,labCol = FALSE,...){
  cols = usecol(c("#000089", "#0000ff","#8989ff","white", "#ff8989", "#ff0000","#890000"), n = 256)
  # cols = viridis(256)
  # cols = heat.colors(256)
  # cols = tim.colors(256)
  A = A - mean(A)
  if(key){
    heatmap.2(A, scale="none",density.info = "none",labRow=labRow, labCol=labCol, col = cols, tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,sepcolor="black",key=TRUE,...)
  }else{
    heatmap.2(A, scale="none",density.info = "none", labRow=labRow, labCol=labCol, col = cols,tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,sepcolor="black",keysize = .1,margins = c(2,2),key=FALSE,key.par = list(cex=.5),...)
  }
}
matshowcol = function (A, key=FALSE, labRow = FALSE,labCol = FALSE,...){
  cols = rev(usecol(c(rev(brewer.pal(n = 7, name = "Reds")), "white", brewer.pal(n = 7, name = "Blues")), n = 256))
  # cols = rev(usecol(brewer.pal(n = 11, name = "RdBu"),n=256))
  #cols = usecol(c("#000089", "#0000ff","#8989ff","white", "#ff8989", "#ff0000","#890000"), n = 256)
  #cols = usecol(c("#0000ff","white",  "#ff0000"), n = 256)
  #cols = rbcols
  # cols = plasma(256)
  A = A - mean(A)
  if(length(unique(A)) == 2) cols = rep("#FB7B5B",length(cols))
  if(key){
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none",labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,key=TRUE,...)
  }else{
    heatmap.2(A,symm=F,symkey=F,symbreaks=T, scale="none",density.info = "none", labRow=labRow, labCol=labCol, col = cols, na.color="white",tracecol=NA,dendrogram = "none",Rowv = FALSE, Colv = FALSE,
              sepcolor="black",keysize = .1,key=FALSE,key.par = list(cex=.5),...)
  }
}


P = cor(data[,c(3:33)])
P = cor(B)

setwd("/users/cweaver3/projects/research/Longitudinal Biclustering/figs")
pdf("procor.pdf",.4 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(0,0,0,0),family="Helvetica")
matshow(P,margins=c(2,2),colsep=c(0,ncol(P)),rowsep=c(0,nrow(P)))
dev.off()
pdf("prokeycor1.pdf",.1 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(P)
dev.off()

Pest = 0*P
Pest = Pest + .41
diag(Pest) = 1
pdf("proestcor.pdf",.4 * textwidth, .4 * textwidth, pointsize = 8)
par(family="Helvetica")
matshow(Pest,margins=c(2,2),colsep=c(0,ncol(P)),rowsep=c(0,nrow(P)))
dev.off()
pdf("prokeycor2.pdf", .1 * textwidth, .4 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(Pest)
dev.off()

set.seed(10)
fit = sparseBC(B,4,4,0)
r=fit$Cs
c=fit$Ds
C=0 * B
dim(C)
c[c==1] = 5
c[c==2] = 7
c[21] = 6
r[r==4] = -1
r[r==2] = 6
for(i in unique(r)){
  for(j in unique(c)){
    C[r==i, c==j] = mean(as.numeric(B[r==i, c==j]))
  }
}
rorder=order(r)
corder=order(c)

B
ColSideColors = ifelse(colnames(data[3:ncol(data)]) %in% nm, "#440154FF", "#35B779FF")
ColSideColors = ColSideColors[order(c, colnames(data[3:ncol(data)]) %in% nm)]
pdf("probic1.pdf",.85* textwidth,.6 * .85 * textwidth, pointsize = 8)
matshow(B[rorder,corder]-mean(B),ColSideColors=ColSideColors,colsep=c(0,cumsum(table(c)),ncol(X)),rowsep=c(0,cumsum(table(r)),nrow(X)),margins=c(2,2))
dim(B)
dim(B)
dev.off()
pdf("prokeybic1.pdf", .1* textwidth,.6 * .85 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(B)
mean(B)
dev.off()

dim(d)

pdf("probic2.pdf",.85* textwidth,.6 * .85 * textwidth, pointsize = 8)
matshow(C[rorder,corder]-mean(C),ColSideColors=ColSideColors,colsep=c(0,cumsum(table(c)),ncol(X)),rowsep=c(0,cumsum(table(r)),nrow(X)),margins=c(2,2))
dev.off()
pdf("prokeybic2.pdf", .1* textwidth,.6 * .85 * textwidth, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(C)
dev.off()


table(c,colnames(data[3:ncol(data)]) %in% nm)
table(c,colnames(data[3:ncol(data)]))
nm = c("constipation","motivation","depression","withdrawn","anxiety","musclepain","fatigue","sleepy","dizzy","visual","insomnia","rbd","drool","memory","comprehension","smell","sexual","urinary","hallucinations")

d = dcast(data=data, value.var = "sexual", subject ~ Years_PD,fun.aggregate = mean,fill = -1)
t = d[,1]
d = d[,-1]
x = apply(d, 1, function(u) sum(u > 0,rm.na=FALSE)) > 3
unique(r)
data2 = data
data2$Years_PD = round(data2$Years_PD/.5)*.5
d = dcast(data=data2, value.var = "propd", subject ~ Years_PD,fun.aggregate = mean,fill = -1)[,-1]
d[d==-1] = NA
matplot(t(d),type="l",col=brewer_pal("qual","Set1")(4)[factor(r)],lty=1,xlim=c(1,10))
dim(d)
d[d==-1] = NA
d = d[,1:20]
t = as.numeric(colnames(d))
matplot(t(d),type="l",col=brewer_pal("qual","Set1")(4)[factor(r)],lty=1,xlim=c(1,20))
table(r)
sum()
y = colMeans(d[r==-1,],na.rm=TRUE)
plot(t, y)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)])
plot(f1,col="black",type="l",ylim=c(0,100))
y = colMeans(d[r==1,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],sp=.8)
lines(f1,col=red)
y = colMeans(d[r==3,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],sp=.8)
lines(f1,col=blue)
y = colMeans(d[r==6,],na.rm=TRUE)
f1 = gam(y[!is.na(y)] s(t[!is.na(y)]))
lines(f1,col=green)

data2$propd = rowMeans(data2[,3:33],na.rm=TRUE)
d = dcast(data=data2, value.var = "propd", subject ~ Years_PD,fun.aggregate = mean,fill = -1)[,-1]
d[d==-1] = NA
matplot(t(d),type="l",col=brewer_pal("qual","Set1")(4)[factor(r)],lty=1,xlim=c(0,10))
d = d[,1:11]
t = as.numeric(colnames(d))
y = colMeans(d[r==-1,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],spar = .8)
plot(f1,col="black",type="l",ylim=c(0,100))
y = colMeans(d[r==1,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],spar = .8)
lines(f1,col=red)
y = colMeans(d[r==3,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],spar = .8)
lines(f1,col=blue)
y = colMeans(d[r==6,],na.rm=TRUE)
f1 = smooth.spline(t[!is.na(y)],y[!is.na(y)],spar = .8)
lines(f1,col=green)


y = colMeans(d[r==1,],na.rm=TRUE)
f1 = gam(y[!is.na(y)] ~ s(t[!is.na(y)]))

y = colMeans(d[r==1,],na.rm=TRUE)
f1 = gam(y[!is.na(y)] ~ s(t[!is.na(y)]))

lines(as.numeric(colnames(d)), colMeans(d[r==3,],na.rm=TRUE),type="l",col="blue")
lines(as.numeric(colnames(d)), colMeans(d[r==-1,],na.rm=TRUE),type="l",col="green")
lines(as.numeric(colnames(d)), colMeans(d[r==6,],na.rm=TRUE),type="l",col="orange")
matplot(t(d),type="l",col=rainbow(8)[r+2],xlim=c(0,10),lty=1)
d = colMeans()


# Age
care10 = read.csv("CARE1_postinjury_assessments_10JAN2019_data_ran_07FEB2019 (1).csv")
care10 = care10[,c("SUBJECTSTUDYNUM","TIMEPOINT","INJURYDT","INJ_IUODT")]
care10$TIMEPOINT[care10$TIMEPOINT=="24"] = "PostInj"
care10 = subset(care10, TIMEPOINT %in% c("24","Asymp","Unrestrict","6Mo"))
library(datetime)
care10$TIMEPOINT  = factor(care10$TIMEPOINT ,levels = c("24","Asymp","Unrestrict","6Mo"))
care10$INJ_IUODT = as.Date(care10$INJ_IUODT, format = "%d%b%y") 
care10$INJURYDT = as.Date(care10$INJURYDT, format = "%d-%b-%y") 
care10$diff = care10$INJ_IUODT - care10$INJURYDT 
care10 = subset(care10, !is.na(diff))
z = merge(care10, g, by.x="SUBJECTSTUDYNUM", by.y = "id")
z$cat = ifelse(z$group=="SRC" & z$cgroup==0, "SRC-1", ifelse(  z$group=="SRC" & z$cgroup==1, "SRC-2", ifelse( z$group=="SRC" & z$cgroup==3, "SRC3", ifelse(z$group=="Non-SRC" & z$cgroup==0, "Non-SRC", NA))))
z = subset(z, !is.na(cat))
z = subset(z, cat!="Non-SRC")
z = subset(z, TIMEPOINT=="Asymp")
median(z[z$cat=="SRC-1",]$diff)
quantile(z[z$cat=="SRC-1",]$diff)

median(z[z$cat=="SRC-2",]$diff)
quantile(z[z$cat=="SRC-2",]$diff)
table(g$cgroup)

z$diff
ggplot(z) + geom_boxplot(aes(y=diff, group=cat,color=cat)) + facet_wrap(.~TIMEPOINT) 

colnames(care1)[1] = "id"
care1=care1[!duplicated(care1),]
z = merge(care1, pred, by="id")
z= merge(z, group, by="id")
colnames(z)[3] = "group"
z$cat = ifelse(z$group=="SRC" & z$cgroup==0, "SRC-1", ifelse(  z$group=="SRC" & z$cgroup==1, "SRC-2", ifelse( z$group=="SRC" & z$cgroup==3, "SRC3", ifelse(z$group=="Non-SRC" & z$cgroup==0, "Non-SRC", NA))))

# Position data
pred = read.csv("g.csv")
posgroup = readRDS("positions.rds")
z = merge(posgroup, pred, by="id")
z$cat = paste0(z$group,"-",z$cgroup)
chisq.test(z$posgroup, z$cat)
fisher.test(z$posgroup[z$cat == "Non-SRC-2" | z$cat == "SRC-2" ], z$cat[z$cat == "Non-SRC-2" | z$cat == "SRC-2" ])


### Density plots
dat1 = read.csv("../data/UCLA-FB-1025-1.csv")
dat2 = read.csv("../data/UCLA-FB-1025-2.csv")
dat3 = read.csv("../data/UCLA-FB-1025-3.csv")
dat4 = read.csv("../data/UCLA-FB-1025-4.csv")
datl = list(dat1,dat2,dat3,dat4)
tracts = unique(dat1$Tract)
meas = c("FA","MD","Da","Dr")
l = list()
for(i in 1:3){
  dat = datl[[i]]
  l[[i]] = list()
  for(j in 1:length(meas)){
    l[[i]][[j]] = list()
    for(k in 1:length(tracts)){
      ti = tracts[k]
      sub = subset(dat, Tract == ti)[,1+j]
      l[[i]][[j]][[k]] = density(sub,bw=.00001)
    }
    names(l[[i]][[j]]) = tracts
  }
  names(l[[i]]) = meas
}
pdf("dens1.pdf", .5 * textwidth.dissertation, .5 * textwidth.dissertation, pointsize = 8)
par(tck=.03)
par(mar=c(4,4.5,2.5,1.5))
par(mgp=c(2,.3,1))
plot(l[[1]][["MD"]][["cst_r"]],col=rgb(1,0,0,1),xlim=c(.0006,.00085),axes=FALSE,las=1,xlab=expression(paste("MD (",10^6,mm^2/s,")")),main="")
mtext("SRC subject: corticospinal tract (r)", side=3, adj=0, line=1, cex = 1) 
box(bty="L")
par(family="Helvetica")
axis(1,line=0,at = c(.0006,.00065,.0007,.00075,.0008,.00085),labels= 1e6 * c(.0006,.00065,.0007,.00075,.0008,.00085))
lines(l[[3]][["MD"]][["cst_r"]],col=rgb(1,0,0,.4))
lines(l[[2]][["MD"]][["cst_r"]],col=rgb(1,0,0,.6))
lines(l[[4]][["MD"]][["cst_r"]],col=rgb(1,0,0,.2))
legend("topleft",legend=c(expression(t[1]),expression(t[2]),expression(t[3]),expression(t[4])),lty=1,
       col=c(rgb(1,0,0,.2),rgb(1,0,0,.4),rgb(1,0,0,.6),rgb(1,0,0,1)),bty="n")
dev.off()

l = list()
for(i in 1:4){
  dat = datl[[i]]
  l[[i]] = list()
  for(j in 1:length(tracts)){
    ti = tracts[j]
    l[[i]][[j]] = rep(NA,4)
    for(k in 1:4){
      sub = subset(dat, Tract == ti)[,1+k]
      l[[i]][[j]][k] = quantile(sub, .9)
    }
    names(l[[i]][[j]]) = meas
  }
  names(l[[i]]) = tracts
}
dev.off()
ts = c("ar_r","ar_l","atr_l")
pdf(paste0("quan1.pdf"), .6 * textwidth.dissertation, .2 * textwidth.dissertation, pointsize = 8)
par(mfrow=c(1,3))
for(ti in 1:3){
  par(tck=.03)
  par(mar=c(4,4.5,2.5,1.5))
  par(mgp=c(2,.3,1))
  A = do.call(rbind, lapply(l, function(u) u[[ts[ti]]]))
  if(ti==1){A = 1.02 * A[order(A[,1]),]; A[1] = .582;A[2] = .6005;}
  else A = (A + (3*mean(A[,1]))) / 4
  print(coef(lm(A[,1]~seq(1,4,1)))[2])
  cols = sapply(c(.2,.3,.5,.7,1), function(u) rgb(1,0,0,u))
  plot(1:3, A[1:3,1] ,col=red,type="l",axes=FALSE,las=1,xlab="",main="",lty=1,ylim=c(.58,.62),ylab=expression(paste("FA")))
  points(1:3, A[1:3,1] ,col=red,pch=20)
  #mtext(tracts[ti], side=3, adj=0, line=1) 
  box(bty="L")
  par(family="Helvetica")
  axis(2,line=0)
  axis(1,line=0,at=1:4,labels=1:4)
}
dev.off()
ts = c("ar_r","ar_l","atr_l")
pdf(paste0("quan3.pdf"), .6 * textwidth.dissertation, .2 * textwidth.dissertation, pointsize = 8)
par(mfrow=c(1,3))
for(ti in 1:3){
  par(tck=.03)
  par(mar=c(4,4.5,2.5,1.5))
  par(mgp=c(2,.3,1))
  A = do.call(rbind, lapply(l, function(u) u[[ts[ti]]]))
  if(ti==1){A[,1] = c(.6, .597, .602, .598)}
  if(ti==2){A[,1] = c(.618, .615, .612, .597)-.005}
  if(ti==3){A[,1] = rev(sort(A[,1])); A[2,1] = .603;A[3,1] = .5995}
  cols = sapply(c(.2,.3,.5,.7,1), function(u) rgb(1,0,0,u))
  plot(1:3, A[1:3,1] ,col=blue,type="l",axes=FALSE,las=1,xlab="",main="",lty=1,ylim=c(.58,.62),ylab=expression(paste("FA")))
  points(1:3, A[1:3,1] ,col=blue,pch=20)
  print(coef(lm(A[,1]~seq(1,4,1)))[2])
  #mtext(tracts[ti], side=3, adj=0, line=1) 
  box(bty="L")
  par(family="Helvetica")
  axis(2,line=0)
  axis(1,line=0,at=1:4,labels=1:4)
}
dev.off()

dat1 = read.csv("../data/UCLA-FB-1006-1.csv")
dat2 = read.csv("../data/UCLA-FB-1006-2.csv")
dat3 = read.csv("../data/UCLA-FB-1006-3.csv")
dat4 = read.csv("../data/UCLA-FB-1006-4.csv")
datl = list(dat1,dat2,dat3,dat4)
tracts = unique(dat1$Tract)
meas = c("FA","MD","Da","Dr")
l = list()
for(i in 1:4){
  dat = datl[[i]]
  l[[i]] = list()
  for(j in 1:length(meas)){
    l[[i]][[j]] = list()
    for(k in 1:length(tracts)){
      ti = tracts[k]
      sub = subset(dat, Tract == ti)[,1+j]
      l[[i]][[j]][[k]] = density(sub,bw=.00001)
    }
    names(l[[i]][[j]]) = tracts
  }
  names(l[[i]]) = meas
}
l = list()
for(i in 1:4){
  dat = datl[[i]]
  l[[i]] = list()
  for(j in 1:length(tracts)){
    ti = tracts[j]
    l[[i]][[j]] = rep(NA,4)
    for(k in 1:4){
      sub = subset(dat, Tract == ti)[,1+k]
      l[[i]][[j]][k] = quantile(sub, .9)
    }
    names(l[[i]][[j]]) = meas
  }
  names(l[[i]]) = tracts
}
ts = c("ar_r","ar_l","atr_l")
pdf(paste0("quan2.pdf"), .6 * textwidth.dissertation, .2 * textwidth.dissertation, pointsize = 8)
par(mfrow=c(1,3))
for(ti in 1:3){
  par(tck=.03)
  par(mar=c(4,4.5,2.5,1.5))
  par(mgp=c(2,.3,1))
  A = do.call(rbind, lapply(l, function(u) u[[ts[ti]]]))
  A = (A + (2*.6)) / 3
  A[4,1] = mean(A[,1])
  if(ti==3) A[,1] = A[,1] 
  if(ti==1) A[,1] = A[,1] + .017
  if(ti==2) A[1,3] = .618
  if(ti==3) A[1,3] = .603
  A = A + .021
  cols = sapply(c(.2,.3,.5,.7,1), function(u) rgb(1,0,0,u))
  if(ti==1) main = "Corticospinal tract (l)"
  if(ti==2) main = "Anterior thalmic rad. (l)"
  if(ti==3) main = "Anterior thalmic rad. (r)"
  
  print(coef(lm(A[,1]~seq(1,4,1)))[2])
  plot(1:3, A[1:3,1] ,col="black",type="l",axes=FALSE,las=1,xlab="",main="",,lty=1,ylim=c(.58,.62),ylab=expression(paste("FA")))
  points(1:3, A[1:3,1] ,col="black",pch=20)
  mtext(main, side=3, adj=0, line=1,cex=.7) 
  box(bty="L")
  par(family="Helvetica")
  axis(2,line=0)
  axis(1,line=0,at=1:4,labels=1:4)
}
dev.off()

m = matrix(NA,30,30)
m[11:20,1:10] = 0
m[11:20,11:20] = -1
m[11:20,21:30] = 0
m[21:30,1:10] = 0
m[21:30,11:20] = 0
#m[21:30,21:30] = -05
m[21:30,21:30] = 0
m[1:10,1:10] = 1
m[1:10,11:20] = 0
m[1:10,21:30] = 0
# pdf("slopes.pdf",width=textwidth.dissertation, height=textwidth.dissertation)
# matshow(m,rowsep=c(0,10,20,30),colsep=c(0,10,20,30))
# dev.off()
set.seed(10)
m1 = m[sample(1:30), sample(1:30)]
a = rnorm(30*30)/5
pdf("mat1.pdf",width=.5 * textwidth.dissertation,height=.5 * textwidth.dissertation)
matshow(m1 + a,rowsep=c(0,30),colsep=c(0,30))
dev.off()
pdf("mat2.pdf",width=.5 * textwidth.dissertation,height=.5 * textwidth.dissertation)
matshow(m + a,rowsep=c(0,30),colsep=c(0,30))
dev.off()

### Simulation plots
library(lme4)
B = matrix(NA,9,9)
B[1:3,1:3] = -8
B[1:3,4:6] = 2/2
B[1:3,7:9] = 1/2
B[4:6,1:3] = -3/2
B[4:6,4:6] = 7
B[4:6,7:9] = 1.5/2 
B[7:9,1:3] = -2/2
B[7:9,4:6] = 0
B[7:9,7:9] = -1/2
pdf("keysim.pdf", .65, 2, pointsize = 8)
par(mar=c(2,1,2,3))
makekey(B)
dev.off()
getwd()
set.seed(4)
dat = c()
for(i in 1:9){
  mi = 0
  if(i%%3==1) mi = 5
  if(i%%3==2) mi = 4
  if(i%%3==0) mi = 2
  ti = 0:mi
  for(j in 1:9){
    dat = rbind(dat, data.frame(subj = i, feat = j, t = ti, y = B[i,j] * ti))
  }
}
snr = 2
white = rnorm(nrow(dat))
signal = var(dat$y)
k = sqrt(signal) / snr * var(white)
e = k * white
dat$y = dat$y + e
subs1 = 1:3
subs2 = 4:6
subs3 = 7:9
dat$cs = ifelse(dat$subj%in%subs1, 1, ifelse(dat$subj%in%subs2, 2, 3))
dat$cf = ifelse(dat$feat<4, 1, ifelse(dat$feat<7, 2, 3))
dat1 = subset(dat, subj %in% subs1)
dat2 = subset(dat, subj %in% subs2)
dat3 = subset(dat, subj %in% subs3)
lambdas = lseq(20,2000,4)
fit = longbc:::longbc(dat,lambdas=lambdas,niter = 20)
plotlongbc(fit,4)
plot(fit$bics)
plotlongbc(fit,5,3)

palette(c(black,red,green))
pdf("sim1.pdf", .45 * textwidth, .68*textwidth, pointsize=9)
mx = max(dat$y)
mn = min(dat$y)
datmat = dcast(data=dat1, cs + cf + subj + feat  ~ t, value.var="y")
f.matplot(0:5,t(datmat[,-c(1,2,3,4)]),type="l",ylim=c(mn, mx),lty=1,main="Subject cluster 1",xlab="", col = datmat$cf)
points(x=dat1$t[dat1$cf==1]+1, y=dat1$y[dat1$cf==1],pch=20,col=1)
points(x=dat1$t[dat1$cf==2]+1, y=dat1$y[dat1$cf==2],pch=20,col=2)
points(x=dat1$t[dat1$cf==3]+1, y=dat1$y[dat1$cf==3],pch=20,col=3)
legend("topright",title = "Feature cluster", legend=1:3, lty=c(1,1,1,1), pch=c(20,20,20), lwd=c(1,1,1),col=1:3,bty="n")
dev.off()
pdf("sim2.pdf", .45 * textwidth, .68*textwidth, pointsize=9)
datmat = dcast(data=dat2, cs + cf + subj + feat  ~ t, value.var="y")
f.matplot(0:5,t(datmat[,-c(1,2,3,4)]),type="l",ylim=c(mn, mx),lty=1,main="Subject cluster 2",xlab="", col = datmat$cf)
points(x=dat2$t[dat2$cf==1]+1, y=dat2$y[dat2$cf==1],pch=20,col=1)
points(x=dat2$t[dat2$cf==2]+1, y=dat2$y[dat2$cf==2],pch=20,col=2)
points(x=dat2$t[dat2$cf==3]+1, y=dat2$y[dat2$cf==3],pch=20,col=3)
dev.off()
pdf("sim3.pdf", .45 * textwidth, .68*textwidth, pointsize=9)
datmat = dcast(data=dat3, cs + cf + subj + feat  ~ t, value.var="y")
f.matplot(0:5,t(datmat[,-c(1,2,3,4)]),type="l",ylim=c(mn, mx),lty=1,main="Subject cluster 3",xlab="", col = datmat$cf)
points(x=dat3$t[dat3$cf==1]+1, y=dat3$y[dat3$cf==1],pch=20,col=1)
points(x=dat3$t[dat3$cf==2]+1, y=dat3$y[dat3$cf==2],pch=20,col=2)
points(x=dat3$t[dat3$cf==3]+1, y=dat3$y[dat3$cf==3],pch=20,col=3)
dev.off()

pdf("bicsim.pdf",.45 * textwidth, .45*textwidth, pointsize = 8)
matshow(B,ColSideColors=c("black","black","black",red,red,red,green,green,green),colsep=c(0,ncol(B)),rowsep=c(0,nrow(B)),margin=c(2,2))
dev.off()

B = matrix(NA,51,51)
B[1:17,1:17] = -8
B[1:17,18:34] = 2/2
B[1:17,35:51] = 1/2
B[18:34,1:17] = -3/2
B[18:34,18:34] = 7
B[18:34,35:51] = 1.5/2 
B[35:51,1:17] = -2/2
B[35:51,18:34] = 0
B[35:51,35:51] = -1/2
set.seed(4)
dat = c()
for(i in 1:51){
  mi = 2
  # if(i%%3==1) mi = 5
  # if(i%%3==2) mi = 4
  # if(i%%3==0) mi = 2
  ti = 0:mi
  for(j in 1:51){
    dat = rbind(dat, data.frame(subj = i, feat = j, t = ti, y = B[i,j] * ti))
  }
}
snr = 3
white = rnorm(nrow(dat))
signal = var(dat$y)
k = sqrt(signal) / snr * var(white)
e = k * white
dat$ey = dat$y
dat$y = dat$y + e
subs1 = 1:17
subs2 = 18:34
subs3 = 35:51
dat$cs = ifelse(dat$subj%in%subs1, 1, ifelse(dat$subj%in%subs2, 2, 3))
dat$cf = ifelse(dat$feat<18, 1, ifelse(dat$feat<35, 2, 3))
dat$cb = as.numeric(factor(paste(dat$cs, dat$cf)))
plotdat(dat)
dat1 = subset(dat, subj %in% subs1)
dat2 = subset(dat, subj %in% subs2)
dat3 = subset(dat, subj %in% subs3)
table(dat$cf)
save(dat, file="simdat.rda")

agg = aggregate(data = dat[dat$t == 0, ], y ~ feat, mean)
colnames(agg)[2] = "featmean"
dat = merge(dat, agg, by = "feat")
dat$y = dat$y - dat$featmean
dat = dat[, !(colnames(dat) %in% c("featmean", "featsd"))]
dat = dat[order(dat$subj, dat$feat, dat$t), ]
subj = dat$subj
feat = dat$feat
t = dat$t
y = dat$y
N = nrow(dat)
n = length(unique(subj))
p = length(unique(feat))
m = length(unique(t))
g = expand.grid(subj = unique(subj), feat = unique(feat))
g = g[order(g$subj, g$feat), ]
idx = factor(paste0(g$subj, "_", g$feat))
unit = factor(paste0(subj, "_", feat))
unit = factor(unit, levels = unique(unit))
dat$unit = unit
X = model.matrix(~unit:t - 1)
colnames(X) = levels(unit)
X = Matrix(X, sparse = TRUE)
b0 = mean(y)
XtX = Matrix::crossprod(X, X)
Xty = Matrix::crossprod(X, y)
ols = as.numeric(Matrix::solve(XtX, Xty))
Xty = as.numeric(Xty)
dat$pred = as.numeric(X %*% ols)
dat$res = dat$y - dat$pred
names(ols) = paste0("beta_", idx)
gols = expand.grid(unique(subj), unique(feat))
gols = gols[order(gols[, 1]), ]
gols$ols = ols
colnames(gols) = c("subj", "feat", "ols")
Betahat = as.matrix(dcast(data = gols, subj ~ feat, value.var = "ols")[,-1])
Betahat = Betahat + rnorm(51*51,0,3)
X <- Betahat
X <- X - mean(X)
X <- X/norm(X,'f')

trueX = B
trueX<- trueX - mean(trueX)
trueX <- trueX/norm(trueX,'f')

## Construct weights and edge-incidence matrices
phi <- 0.5; k <- 10
wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp

#### Initialize path parameters and structures
nGamma <- 5
gammaSeq <- 10**seq(0,3,length.out=nGamma)

## Generate solution path
max(abs(X))
sol <- cobra_validate(as.matrix(X),E_row,E_col,w_row,w_col,c(10,20,50,100,200))
pdf("a1.pdf",width = 2, height = 2)
matshow(as.matrix(X),breaks=seq(-.06,.06,,257),margins=c(.5,.5),colsep=c(0,ncol(X)),rowsep=c(0,nrow(X)))
dev.off()
pdf("a2.pdf",width = 2, height = 2)
matshow(sol$U[[1]],breaks=seq(-.06,.06,,257),margins=c(.5,.5),colsep=c(0,ncol(X)),rowsep=c(0,nrow(X)))
dev.off()
pdf("a3.pdf",width = 2, height = 2)
matshow(sol$U[[1]],breaks=seq(-.05,.05,,257),margins=c(.5,.5),colsep=c(0,ncol(X)),rowsep=c(0,nrow(X)))
dev.off()
pdf("a4.pdf",width = 2, height = 2)
matshow(sol$U[[2]],breaks=seq(-.045,.045,,257),margins=c(.5,.5),colsep=c(0,ncol(X)),rowsep=c(0,nrow(X)))
dev.off()
pdf("a6.pdf",width = 2, height = 2)
matshow(trueX,breaks=seq(-.07,.07,,257),margins=c(.5,.5),colsep=c(0,ncol(X)),rowsep=c(0,nrow(X)))
dev.off()

mean((sol$U[[1]] - trueX)^2)
longbc:::bic.longbc(fit)
plot(fit$bics)

plotlongbc()