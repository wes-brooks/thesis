library(lagr)
library(mgcv)
library(dplyr)
library(doMC)

registerDoMC(7)

set.seed(11181982)

#Set constants:
B = 20
n = 100

#Simulate data:
f0 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f1 = function(x) {cos(x)}
f2 = function(x) {0.1*sin(x)}
tt = seq(0, 5, len=n)


f0.second.derivative = function(x) {-2.4*x^2 + 6*x + 1.4}
f1.second.derivative = function(x) {-cos(x)}
f2.second.derivative = function(x) {-sin(x)}

fitted = list()
resid = list()
coefs = list()

empirical.bias.0 = list()
empirical.bias.1 = list()
empirical.bias.2 = list()

for (b in 1:B) {
print(b)
    X1 = rnorm(n, mean=3, sd=2)
    X2 = rnorm(n, mean=3, sd=2)
    y = f0(tt) + X1*f1(tt) + X2*f2(tt) + rnorm(n)
    df = data.frame(y, x1=X1, x2=X2, t=tt)

    m = lagr(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAIC', kernel=epanechnikov, bw.type='knn', bw=0.2, verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    
    W = list()
    for (i in 1:n) {
        h = m[['fits']][[i]][['bw']]
        w = epanechnikov(abs(tt-tt[i]), h)
        W[[i]] = diag(w)
    }
    
    fitted[[b]] = sapply(m[['fits']], function(x) tail(x[['fitted']],1))
    resid[[b]] = df$y - fitted[[b]]
    coefs[[b]] = t(sapply(m[['fits']], function(x) x[['model']][['beta']][,ncol(x[['model']][['beta']])]))

    empirical.bias.0[[b]] = c(0, diff(coefs[[b]][,4])) * n / 5 * sapply(m[['fits']], function(x) x[['bw']]) / 10
    empirical.bias.1[[b]] = c(0, diff(coefs[[b]][,5])) * n / 5 * sapply(m[['fits']], function(x) x[['bw']]) / 10
    empirical.bias.2[[b]] = c(0, diff(coefs[[b]][,6])) * n / 5 * sapply(m[['fits']], function(x) x[['bw']]) / 10
}
bias.0 = f0.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]) / 10
bias.1 = f1.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]) / 10
bias.2 = f2.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]) / 10

empirical.bias = list()
for (b in 1:B) {
    empirical.bias[[b]] = cbind(empirical.bias.0[[b]], empirical.bias.1[[b]], empirical.bias.2[[b]])
}

sapply(empirical.bias.0, c) %>% rowMeans -> eb0


#sqrt(n*sapply(m[['fits']], function(x) x[['bw']]) / 5) -> bias.factor

#Plot multiple realizations of the means: resampling should look something like this(?)
yy = range(c(f2(tt), sapply(coefs, function(x) x[,3])))
plot(x=tt,y=f2(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
for (b in 1:B) {
    par(new=TRUE)
    #Subtract the bias as it appears from the rate of change of the first derivative:
    plot(x=tt, y=coefs[[b]][,3], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
}
par(new=TRUE)
plot(x=tt, y=f2(tt) + bias.2, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)




#Linked bootstrap draws:
beta.star.1 = list()

for (j in 1:B) {
    cat(paste("j is ", j, "\n", sep=""))
    beta.star.1[[j]] = matrix(NA,6,0)
    y.hat.seed = as.matrix(rnorm(6))
    
    for (i in 1:n) {
        Sigma = with(summary(m[['fits']][[i]][['model']][['adamodel']]), cov.unscaled * dispersion)[1:6,1:6]
        beta.star.1[[j]] = cbind(beta.star.1[[j]], (coefs[[B]][i,] + (sqrtm(Sigma) %*% y.hat.seed)))
    }
    beta.star.1[[j]] = beta.star.1[[j]][1:3,] - t(empirical.bias[[B]])
}


#Plot draws from the linked bootstrap: TOO SMOOTH
yy2 = sapply(beta.star.1, function(x) range(x[2,])) %>% range
plot(x=tt,y=f1(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy2, type='l')
par(new=TRUE)
plot(x=tt, y=coefs[[B]][,2], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue')
par(new=TRUE)
plot(x=tt, y=f1(tt) + bias.1, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue', lty=2)
for (b in 1:B) {
    par(new=TRUE)
    plot(x=tt, y=beta.star.1[[b]][2,], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='red')
}



#Now a full distribution bootstrap (too rough?):
Sigma = list()
for (i in 1:n) {
    Sigma[[i]] = with(summary(m[['fits']][[i]][['model']][['adamodel']]), cov.unscaled * dispersion)[1:3,1:3]
}

SS = bdiag(Sigma)
for (i in 1:(n-1)) {
    X.i = cbind(1, X1, X2)
    
    for (j in (i+1):n) {
        X.j = cbind(1, X1, X2)
        scale = sqrt(summary(m[['fits']][[i]][['model']][['adamodel']])$dispersion *
                         summary(m[['fits']][[j]][['model']][['adamodel']])$dispersion)
        
        unscaled = summary(m[['fits']][[i]][['model']][['adamodel']])$cov.unscaled[1:3,1:3] %*%
            t(X.i) %*% sqrt(W[[i]] %*% W[[j]]) %*% 
            X.j %*% summary(m[['fits']][[j]][['model']][['adamodel']])$cov.unscaled[1:3,1:3]
        
        SS[((i-1)*3+1):(3*i), ((j-1)*3+1):(3*j)] =
            SS[((j-1)*3+1):(3*j), ((i-1)*3+1):(3*i)] = scale * unscaled
    }
}

SS.posdef = SS %*% t(SS)
SS.eigen = eigen(SS.posdef)

sqSS = SS.eigen$vectors %*% diag(SS.eigen$values^(1/4)) %*% t(SS.eigen$vectors)

#Draw from a joint parametric distribution:
beta.star.3 = list()
beta.star.4 = list()
for (k in 1:B) {
    #beta.star.3[[k]] = (t(coefs) + matrix(sqSS %*% rep(rnorm(6),n), 6, n))
    beta.star.4[[k]] = (t(coefs[[B]][,1:3]) + matrix(sqSS %*% rnorm(3*n), 3, n))
}

yy2 = sapply(beta.star.4, function(x) range(x[1,])) %>% range
plot(x=tt,y=f0(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy2, type='l')
par(new=TRUE)
plot(x=tt, y=coefs[[B]][,1], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue')
for (b in 1:B) {
    par(new=TRUE)
    plot(x=tt, y=beta.star.4[[b]][1,], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='red')
}



g1 = gam(resid[[B]]~s(tt, k=50))
s2 = g1$sig2

plot(x=tt, y=resid[[B]], bty='n', xlab='t', ylab='residual')
par(new=TRUE)
plot(x=tt, y=fitted(g1), type='l', ylim=range(resid), ann=FALSE, xaxt='n', yaxt='n', bty='n')
