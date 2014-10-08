library(mvtnorm)
library(dplyr)
library(glmnet)

sigma = matrix(0.05,5,5)
diag(sigma) = rep(1,5)
x = rmvnorm(n=1000, mean=c(0,0,0,0,0), sigma=sigma)
y = rnorm(1000, mean=x[,1] + 0.3*x[,2], sd=1)
p = ncol(x)

beta.hat = lm(y~x)$coef[2:(p+1)]
x.scaled = sweep(x, 2, beta.hat, '/') %>% as.matrix

m = glmnet(y=y, x=x.scaled)
cc = as.matrix(coef(m))
fitted = cbind(1,x.scaled) %*% cc

AIC = apply(fitted, 2, function(x) sum((x-y)**2)) + 2*m$df
k = which.min(AIC)

fitted = fitted[,k]
beta = cc[,k]/c(1,beta.hat)
r = y - fitted

for (i in 1:100000) {
    rr = sample(r, replace=TRUE)
    yy = fitted + rr - mean(rr)
    
    beta.hat.boot = lm(yy~x)$coef[2:(p+1)]
    x.boot = sweep(x, 2, beta.hat.boot, '/') %>% as.matrix
    m.boot = glmnet(y=yy, x=x.boot)
    
    fit.boot = cbind(1,x.boot) %*% as.matrix(coef(m.boot))
    AIC.boot = apply(fit.boot, 2, function(x) sum((x-y)**2)) + 2*m.boot$df
    k.boot = which.min(AIC.boot)
    
    beta = cbind(beta, coef(m.boot)[,k.boot]/c(1,beta.hat.boot))
}

for (j in 1:6) {
    ee = ecdf(beta[j,])
    
}
