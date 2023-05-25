load('GSE72680_422.Rdata')

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

ex <- expr
ex_m <- log2(ex/(1-ex))   ##log-transformed methylation M-values 
ex_m=ex_m[,-ID.m] # remove missing values in y and x

anno=read.csv("annotation.csv", header = TRUE, row.names = 1, sep = ",", quote = "\"",
           dec = ".", fill = TRUE, comment.char = "")
row.names1 <- rownames(ex_m)
row.names2 <- rownames(anno)
matched_indices <- match(row.names2, row.names1)
ex_m <- ex_m[matched_indices, ] # methylation in gene promoter region
#average 
ex_m=t(ex_m)

fit0 <- lm(y ~ x+cov[,1:3]+as.factor(cov[,4]), na.action = na.exclude) # total effect
summary(fit0)

fit1 = lm(ex_m[,1]~x+cov[,1:3]+as.factor(cov[,4])+cellp, na.action = na.exclude) # mediation model
summary(fit1)

fit2 = lm(y~ex_m[,1]+x+cov[,1:3]+as.factor(cov[,4])+cellp, na.action = na.exclude) # outcome model
summary(fit2)

tme=fit1$coefficients[2]*fit2$coefficients[2] # total mediation effect a*b
tme 

