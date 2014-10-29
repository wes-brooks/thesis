library(lagr)

set.seed(12345)
source("scratch/inference/simulate-data.r")

source("scratch/inference/read-trace.r")

#Set the location to fit the model:
fl = sim[106,c('loc.x','loc.y')]

coefs = matrix(NA,0,5)

for (h in seq(0.18, 0.25, len=15)) {
for (i in 1:10) {
print(i)

#Generate the bootstrap draws:
indx = sample(1:(N^2), replace=TRUE)


m = lagr(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), fit.loc=fl, longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005, bootstrap.index=indx)

beta = m[['model']][[1]][['model']][['beta']][1:5,]
models = as.factor(apply(beta, 2, function(x) paste(as.integer(x!=0), collapse="")))
aic.trace = m[['model']][[1]][['model']][['results']][['AIC']]

best = vector()
for (l in levels(models)) {
    best = c(best, which.min(ifelse(models==l, aic.trace, Inf)))
}

coefs = rbind(coefs, apply(beta[,best], 1, function(x) weighted.mean(x, w=exp(-2*aic.trace[best]))))
}
}