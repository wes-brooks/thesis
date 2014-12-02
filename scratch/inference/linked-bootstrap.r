library(lagr)
library(MASS)
library(expm)
library(mgcv)
library(Matrix)

B = 1
n = 50
h = 1.5

f1 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f2 = function(x) {cos(x)}
tt = seq(0, 5, len=n)
X = rnorm(n, mean=3, sd=2)
y = f1(tt) + X*f2(tt) + rnorm(n)

m = list()
fitted = vector()
resid = vector()
estimate = vector()
W = list()
fitted.coefs = matrix(NA, 0, 4)
for (i in 1:n) {
    w = epanechnikov(abs(tt-tt[i]), h)
    W[[i]] = diag(w)
    X.loc = cbind(tt-tt[i], X, (tt-tt[i])*X)
    m[[i]] = lm(y~X.loc, weights=w)
    fitted = c(fitted, m[[i]]$fitted[i])
    resid = c(resid, m[[i]]$resid[i])
    estimate = c(estimate, m[[i]]$coef[3])
    
    fitted.coefs = rbind(fitted.coefs, coef(m[[i]]))
}

#Draws:
y.hat1 = list()
y.hat2 = list()
for (j in 1:B) {
    cat(paste("j is ", j, "\n", sep=""))
    y.hat1[[j]] = vector()
    y.hat.seed = as.matrix(rnorm(4))
    y.hat2[[j]] = vector()
    for (i in 1:n) {
        Sigma = with(summary(m[[i]]), cov.unscaled * sigma^2)
        y.hat1[[j]] = c(y.hat1[[j]], mvrnorm(1, mu=m[[i]]$coef, Sigma=Sigma)[1])
        y.hat2[[j]] = c(y.hat2[[j]], (m[[i]]$coef + (sqrtm(Sigma) %*% y.hat.seed))[1])
    }
}


# i = 30
# for(i in 1:40) {
#     #yy = range(c(y.hat1[[i]], y.hat2[[i]], f2(tt)))
#     #plot(y.hat1[[i]], type='l', col='red', ylim=yy)
#     par(new=TRUE)
#     plot(y.hat2[[i]], type='l', col='blue', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
# }
# par(new=TRUE)
# plot(m[[1]]$fitted, type='l', col='red', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
# 
# 
# #Smooth the residuals with a GAM:
# g1 = gam(resid~s(tt, k=50))
# s2 = g1$sig2


#Draw coefficients two ways:
#1) independent bootstrap
#2) Linked bootstrap
beta.hat1 = list()
beta.hat2 = list()
for (j in 1:B) {
    cat(paste("j is ", j, "\n", sep=""))
    beta.hat1[[j]] = matrix(NA,2,0)
    y.hat.seed = as.matrix(rnorm(4))
    beta.hat2[[j]] = matrix(NA,2,0)
    for (i in 1:n) {
        Sigma = with(summary(m[[i]]), cov.unscaled * sigma^2)
        beta.hat1[[j]] = cbind(beta.hat1[[j]], mvrnorm(1, mu=m[[i]]$coef, Sigma=Sigma)[1:2])
        beta.hat2[[j]] = cbind(beta.hat2[[j]], (m[[i]]$coef + (sqrtm(Sigma) %*% y.hat.seed))[1:2])
    }
}


Sigma = list()
for (i in 1:n) {
    Sigma[[i]] = with(summary(m[[i]]), cov.unscaled * sigma^2)
}

SS = bdiag(Sigma)
for (i in 1:(n-1)) {
    X.i = cbind(1, tt-tt[i], X, (tt-tt[i])*X)
    
    for (j in (i+1):n) {
        X.j = cbind(1, tt-tt[j], X, (tt-tt[j])*X)
        scale = summary(m[[i]])$sigma * summary(m[[j]])$sigma
        
        unscaled = summary(m[[i]])$cov.unscaled %*% t(X.i) %*% W[[i]] %*% W[[j]] %*% X.j %*% summary(m[[j]])$cov.unscaled
        SS[((i-1)*4+1):(4*i), ((j-1)*4+1):(4*j)] = SS[((j-1)*4+1):(4*j), ((i-1)*4+1):(4*i)] = scale * unscaled
    }
}

sqSS = sqrtm(SS)
coef.draw = sqSS %*% rnorm(4*n)

#yy = range(c(beta.hat1[[i]][2,], beta.hat2[[i]][2,], f2(tt)))
#plot(f2(tt), type='l', col='black', ylim=yy)
#par(new=TRUE)
#plot(beta.hat1[[i]][2,], type='l', col='red', lty=2, ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
#par(new=TRUE)
#plot(beta.hat2[[i]][2,], type='l', col='blue', lty=3, ylim=yy, ann=FALSE, xaxt='n', yaxt='n')

#for(i in 1:40) {
#    #yy = range(c(beta.hat1[[i]][2,], beta.hat2[[i]][2,], f2(tt)))
#    #plot(beta.hat1[[i]][2,], type='l', col='red', ylim=yy)
#    #i=9
#    par(new=TRUE)
#    plot(beta.hat2[[i]][2,], type='l', col='blue', ylim=yy, ann=FALSE, xaxt='n', yaxt='n')
#}
#par(new=TRUE)
#plot(f2(tt), type='l', col='black', lwd=2, ylim=yy, ann=FALSE, xaxt='n', yaxt='n')



