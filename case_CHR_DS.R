################### functions ##################
f = function(u,v){
  return(u + v)
}

fdr = function(t, M){
  up = sum(M<(-1*t))
  down = max(1, sum(M>t))
  return(up/down)
}

################### Args ###################
require("getopt", quietly=TRUE)

spec = matrix(c(
    "Chr", "c", 1, "integer"
), byrow=TRUE, ncol=4)

opt = getopt(spec);
cat('Chr:', opt$Chr,'.\n')


###################### data preparation ####################
load('GSE72680_422.Rdata')
#y: response
#x: treatment/exposure
#cellp, cov: pre-treatment covariates
#ex_m: mediations
y <- pdata$PSS #PTSD_Dep
x <- pdata$sexual.abuse;
x <- as.numeric(as.character(x));
x[which(x==0)]=1
x=x-1;

#x=pdata$child.trauma.score #this is the exposure or treatment
ID.m=which(is.na(y) | is.na(x))
pdata2=pdata[-ID.m,] # remove missing values in y and x
y=y[-ID.m]
x=x[-ID.m]

cov <- pdata2[6:9]# race is a nominal variable, treated as a factor
cov <- do.call(rbind, cov)
cov=t(cov)

cellp=pdata2[,-c(1:10)] # these are cell compositions need to be adjusted as covariates whenever methylation is included in the model.   
cellp <- do.call(rbind, cellp)
cellp=t(cellp)

## preparation of mediations
ex <- expr
ex_m <- log2(ex/(1-ex))   ##log-transformed methylation M-values 
ex_m=ex_m[,-ID.m] # remove missing values in y and x

anno=read.csv("annotation.csv", header = TRUE, row.names = 1, sep = ",", quote = "\"",
              dec = ".", fill = TRUE, comment.char = "")
chr = rownames(anno)[which(anno$CHR==opt$Chr)] # mediations in CHR
ex_m_chr = ex_m[which(rownames(ex_m) %in% chr),]
p = dim(ex_m_chr)[1] # number of mediators
cat("Number of mediators:", p, '.\n')
chr_name = rownames(ex_m_chr)
saveRDS(chr_name, file = paste('chr_names/name_CHR',as.character(opt$Chr), '.rds', sep=''))
cat('Chr names saved. \n')

## final removal of NA
y = y[!is.na(cov[,3])]
x = x[!is.na(cov[,3])]
cellp = cellp[!is.na(cov[,3]),]
ex_m_chr = ex_m_chr[,!is.na(cov[,3])]
cov = cov[!is.na(cov[,3]),]

## covert every variable to matrix
library(fastDummies)
y.m = matrix(y, ncol=1)
t.m = matrix(x, ncol=1)
cellp.m = as.matrix(cellp)
cov.dum = dummy_cols(cov, select_columns = c('sex', 'race'))[, c(1,3,5,6,7,8,9,10)]
cov.m = as.matrix(cov.dum)
cov.m = cbind(cellp.m, cov.m)
med.m = matrix(ex_m_chr, ncol=ncol(ex_m_chr))
med.m = t(med.m)


################### DS #############################
library(MASS)
library(stats)
library(glmnet)

S1.list = list()

## split data into equal parts
n = length(y.m)
ids = c(1:n)
for(rep in 1:100){
D1 = sample(ids, size = n/2)
D2 = setdiff(ids, D1)

y1 = y.m[D1]
y2 = y.m[D2]

t1 = t.m[D1]
t2 = t.m[D2]

inter.1 = matrix(rep(1,n/2), ncol = 1) #intercept
inter.full = matrix(rep(1,n), ncol = 1)

cov1 = cov.m[D1,]
cov2 = cov.m[D2,]

m1 = med.m[D1,]
m2 = med.m[D2,]

## phase1: Lasso on outsome model Y
X.full1 = cbind(inter.1, t1, cov1, m1)
p.fac1 = rep(1, dim(X.full1)[2])
p.fac1[1:15] = 0

cv.lasso.1 = cv.glmnet(X.full1, y1, penalty.factor = p.fac1)
coef.lasso1.min = coef(cv.lasso.1$glmnet.fit, s = cv.lasso.1$lambda.min, exact = F)

ind.lasso.raw = which(coef.lasso1.min!=0)
# estimation of beta1 we need
beta1.e1 = coef.lasso1.min[ind.lasso.raw[ind.lasso.raw>16]]

ind.lasso = ind.lasso.raw[ind.lasso.raw>16]
ind.lasso = ind.lasso - 16#index of active mediators, p0 mediators are selected
p1.hat = length(ind.lasso)

if(p1.hat != 0){
# Phase2: Do regression on D2 and ind.lasso for gamma and beta
# Mediator model
M.S1 = med.m[, ind.lasso]
model.m = lm(M.S1~t.m+cellp.m+as.factor(cov[,2])+as.factor(cov[,4])+cov[,c(1,3)])
if(p1.hat == 1){
  gamma.e1 = model.m$coefficients[2]
}else{
  gamma.e1 = model.m$coefficients[2,]
}

# Outcome model
M2.S1 =  m2[, ind.lasso]
cov.2 = cov[D2,]
cellp.2 = cellp[D2,]
model.o = lm(y2~M2.S1+t2+cellp.2+as.factor(cov.2[,2])+as.factor(cov.2[,4])+cov.2[,c(1,3)])
beta2.e1 = model.o$coefficients[(2:(1+p1.hat))]

# generate mirror statistics
min.1 = c()
min.2 = c()
for(i in 1:length(gamma.e1)){
  min.1 = c(min.1, sign(gamma.e1[i] * beta1.e1[i]) * min(abs(gamma.e1[i]), abs(beta1.e1[i])))
  min.2 = c(min.2, sign(gamma.e1[i] * beta2.e1[i]) * min(abs(gamma.e1[i]), abs(beta2.e1[i])))
}

Mirror = c()
for(i in 1:length(min.1)){
  Mirror = c(Mirror, sign(min.1[i]*min.2[i])*f(abs(min.1[i]), abs(min.2[i])))
}

Mirror.abs = abs(Mirror)
fdr.M = c()
for(i in Mirror.abs){
  fdr.M = c(fdr.M, fdr(i, Mirror))
}
cutoffs = Mirror.abs[fdr.M<=0.05]
cutoff = min(cutoffs)

S1.final = ind.lasso[which(Mirror > cutoff)]#final set of selected values.
S0.final = setdiff(1:p, S1.final)

if(length(S1.final)!=0){
  S1.list[[rep]] = S1.final
}else{
  cat('0 mediator is selected by DS.\n')
}
}else{
  cat('0 mediator is selected by LASSO.\n')
}
}

S1.list = Filter(Negate(is.null), S1.list)
saveRDS(S1.list, file = paste('selection/CHR',as.character(opt$Chr), '_DS.rds', sep=''))