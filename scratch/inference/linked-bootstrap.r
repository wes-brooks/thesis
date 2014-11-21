f1 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f2 = function(x) {cos(x)}
tt = seq(0, 5, len=100)
X = rnorm(100, mean=3, sd=2)

y = f1(tt) + X*f2(tt) + rnorm(100)
h = 1.5

m = list()
fitted = vector()
resid = vector()
estimate = vector()
for (i in 1:100) {
    w = epanechnikov(abs(tt-tt[i]), h)
    X.loc = cbind(X, tt-tt[i], (tt-tt[i])*X)
    m[[i]] = lm(y~X.loc, weights=w)
    fitted = c(fitted, m[[i]]$fitted[i])
    resid = c(resid, m[[i]]$resid[i])
    estimate = c(estimate, m[[i]]$coef[3])
}

#Draws:
y.hat1 = list()
y.hat2 = list()
for (j in 1:40) {
    cat(paste("j is ", j, "\n", sep=""))
    y.hat1[[j]] = vector()
    y.hat.seed = as.matrix(rnorm(2))
    y.hat2[[j]] = vector()
    for (i in 1:100) {
        Sigma = with(summary(m[[i]]), cov.unscaled * sigma^2)
        y.hat1[[j]] = c(y.hat1[[j]], mvrnorm(1, mu=m[[i]]$coef, Sigma=Sigma)[1])
        y.hat2[[j]] = c(y.hat2[[j]], (m[[i]]$coef + (sqrtm(Sigma) %*% y.hat.seed))[1])
    }
}


i = 30
for(i in 1:40) {
    #yy = range(c(y.hat1[[i]], y.hat2[[i]], f(tt)))
    #plot(y.hat1[[i]], type='l', col='red', ylim=yy)
    par(new=TRUE)
    plot(y.hat2[[i]], type='l', col='blue', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
}
par(new=TRUE)
plot(yhat$fitted, type='l', col='red', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')


#Smooth the residuals with a GAM:
g1 = gam(resid~s(tt, k=50))
s2 = g1$sig2


#Draw coefficients two ways:
#1) independent bootstrap
#2) Linked bootstrap
beta.hat1 = list()
beta.hat2 = list()
for (j in 1:40) {
    cat(paste("j is ", j, "\n", sep=""))
    beta.hat1[[j]] = matrix(NA,2,0)
    y.hat.seed = as.matrix(rnorm(4))
    beta.hat2[[j]] = matrix(NA,2,0)
    for (i in 1:100) {
        Sigma = with(summary(m[[i]]), cov.unscaled * sigma^2)
        beta.hat1[[j]] = cbind(beta.hat1[[j]], mvrnorm(1, mu=m[[i]]$coef, Sigma=Sigma)[1:2])
        beta.hat2[[j]] = cbind(beta.hat2[[j]], (m[[i]]$coef + (sqrtm(Sigma) %*% y.hat.seed))[1:2])
    }
}



i = 30
for(i in 1:40) {
    #yy = range(c(beta.hat1[[i]][2,], beta.hat2[[i]][2,], f2(tt)))
    #plot(beta.hat1[[i]][2,], type='l', col='red', ylim=yy)
    i=9
    par(new=TRUE)
    plot(beta.hat2[[i]][2,], type='l', col='blue', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
}
par(new=TRUE)
plot(yhat$fitted, type='l', col='red', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')



