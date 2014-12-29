## @knitr simulate-data
#Set constants:
B1 = 20
n = 100
tt = seq(0, 5, len=n)

#Coefficient functions:
f0 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f1 = function(x) {0.2*x^2 - x + cos(x)}
f2 = function(x) {(x+1)^(-1) - 1/6}

f0.second.derivative = function(x) {-2.4*x^2 + 6*x + 1.4}
f1.second.derivative = function(x) {0.4 - cos(x)}
f2.second.derivative = function(x) {2*(x+1)^(-3)}

f0.fourth.derivative = function(x) {-4.8}
f1.fourth.derivative = function(x) {cos(x)}
f2.fourth.derivative = function(x) {24*(x+1)^(-5)}

#Create objects to store results:
fitted = list()
resid = list()
coefs = list()
bw = list()

#Generate the raw data sets (not bootstrap)
set.seed(123)
for (b in 1:B1) {
    #simulate data
    X1 = rnorm(n, mean=3, sd=2)
    X2 = rnorm(n, mean=3, sd=2)
    y = f0(tt) + X1*f1(tt) + X2*f2(tt) + rnorm(n)
    df = data.frame(y, x1=X1, x2=X2, t=tt)
    
    bw[[b]] = lagr.tune(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc',
                   kernel=epanechnikov, bw.type='dist', bwselect.method='AIC', tol.bw=0.1, verbose=FALSE,
                   lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    
    #Fit a VCR model to the simulated data by LAGR
    m = lagr(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc', bw=bw, verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    
    #Get the observation weights for each local fit in the VCR model:
    W = list()
    for (i in 1:n) {
        h = m[['fits']][[i]][['bw']]
        w = epanechnikov(abs(tt-tt[i]), h)
        W[[i]] = diag(w)
    }
    
    #Store the fitted values, residuals, and coefficient estimates:
    fitted[[b]] = sapply(m[['fits']], function(x) tail(x[['fitted']],1))
    resid[[b]] = df$y - fitted[[b]]
    coefs[[b]] = t(sapply(m[['fits']], function(x) x[['model']][['beta']][,ncol(x[['model']][['beta']])]))
}


## @knitr asymptotic-bias

#Compute asymptotic bias from the true coefficient functions:
bias.0 = f0.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.1 = f1.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.2 = f2.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10

#Asymptotic bias of the second order, for parametric bootstrap estimates of the coefficients
bias.0.2 = f0.fourth.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.1.2 = f1.fourth.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.2.2 = f2.fourth.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10






## @knitr empirical-bias

#Empirical bias:
empirical.bias.0 = list()
empirical.bias.1 = list()
empirical.bias.2 = list()
empirical.bias = list()

for (b in 1:B1) {
    #Empirical second derivative of beta.0
    sl = coefs[[b]][,4]
    cc.0 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.0 = rbind(cc.0, m2$coef)
    }
    
    empirical.bias.0[[b]] = cc.0[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
    
    #Empirical second derivative of beta.1
    sl = coefs[[b]][,5]
    cc.1 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.1 = rbind(cc.1, m2$coef)
    }
    
    empirical.bias.1[[b]] = cc.1[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
    
    #Empirical second derivative of beta.2
    sl = coefs[[b]][,6]
    cc.2 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.2 = rbind(cc.2, m2$coef)
    }
    
    empirical.bias.2[[b]] = cc.2[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
    
    empirical.bias[[b]] = cbind(empirical.bias.0[[b]], empirical.bias.1[[b]], empirical.bias.2[[b]])
}





## @knitr plot-realization
layout(matrix(1:3,1,3))

#Plot a realization of the fitted model: \beta_0
yy = range(c(f0(tt), sapply(coefs, function(x) x[,1])))
plot(x=tt,y=f0(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,1])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
par(new=TRUE)
plot(x=tt, y=f0(tt) + bias.0, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,1])) - rowMeans(sapply(empirical.bias, function(b) b[,1])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='green', lty=2, lwd=2)

#Plot a realization of the fitted model: \beta_1
yy = range(c(f1(tt), sapply(coefs, function(x) x[,2])))
plot(x=tt, y=f1(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,2])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
par(new=TRUE)
plot(x=tt, y=f1(tt) + bias.1, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,2])) - rowMeans(sapply(empirical.bias, function(b) b[,2])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='green', lty=2, lwd=2)

#Plot a realization of the fitted model: \beta_2
yy = range(c(f2(tt), sapply(coefs, function(x) x[,3])))
plot(x=tt,y=f2(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,3])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
par(new=TRUE)
plot(x=tt, y=f2(tt) + bias.2, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,3])) - rowMeans(sapply(empirical.bias, function(b) b[,3])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='green', lty=2, lwd=2)




## @knitr linked-bootstrap-sample
B = 10

#Linked bootstrap draws to generate the resampled response:
beta.star.1 = list()
Y.star.1 = list()
X = as.matrix(cbind(1, X1, X2))

for (j in 1:B) {
    cat(paste("j is ", j, "\n", sep=""))
    beta.star.1[[j]] = matrix(NA,6,0)
    y.hat.seed = as.matrix(rnorm(6))
    
    for (i in 1:n) {
        Sigma = with(summary(m[['fits']][[i]][['model']][['adamodel']]), cov.unscaled * dispersion)[1:6,1:6]
        beta.star.1[[j]] = cbind(beta.star.1[[j]], (coefs[[B1]][i,] + (sqrtm(Sigma) %*% y.hat.seed)))
    }
    beta.star.1[[j]] = t(beta.star.1[[j]][1:3,]) - empirical.bias[[B1]]
    Y.star.1[[j]] = sapply(1:n, function(k) t(X[k,]) %*% beta.star.1[[j]][k,]) + rnorm(n, 0, sd=sd(resid[[B1]]))
}






## @knitr linked-bootstrap-estimate

#Run estimation on the parametric bootstrap draws:
coefs.boot = list()
conf.zero = list()
df.b = df
for (b in 1:B) {
    print(b)
    df.b$y = Y.star.1[[b]]
    
    m.b = lagr(y~x1+x2, data=df.b, family='gaussian', coords='t', varselect.method='wAICc', kernel=epanechnikov, bw.type='dist', bw=1, verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    coefs.boot[[b]] = t(sapply(m.b$fits, function(x) x$coef))
    conf.zero[[b]] = t(sapply(m.b$fits, function(x) x$conf.zero))
}


#sd by Efron method (nonparametric delta method):
sd.boot = list()
zero.boot = list()
for (i in 1:n) {    
    t.b = sapply(coefs.boot, function(x) x[i,]) %>% t
    t.mean = matrix(rep(colMeans(t.b),each=B),B,3)
    B.b = (sapply(beta.star.1, function(x) x[i,]) %>% t)[1:B,]
    
    z.b = sapply(conf.zero, function(x) x[i,-1]) %>% t
    z.mean = matrix(rep(colMeans(z.b),each=B),B,2)
    coz.b = t(B.b) %*% (z.b-z.mean) / B
    
    cov.b = t(B.b) %*% (t.b-t.mean) / B
    V.b = t(B.b) %*% B.b / B
    sd.boot[[i]] = sqrtm(t(cov.b) %*% solve(V.b) %*% cov.b)
    
    zero.boot[[b]] = sqrtm(t(coz.b) %*% solve(V.b) %*% coz.b)
}


## @knitr linked-bootstrap-plot-quantiles

#png("/Users/wesley/git/thesis/writeup/paraboot/figures/paraboot-est.png", width=7, height=3, units='in', res=150)
layout(matrix(1:3,1,3))
#confidence from paraboot:
sapply(coefs.boot, function(x) x[,1]) %>% range -> yy
sapply(coefs.boot, function(x) x[,1]) %>% apply(1, function(y) quantile(y, 0.05)) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[0]), ylim=yy, bty='n')
par(new=TRUE)
sapply(coefs.boot, function(x) x[,1]) %>% apply(1, function(y) quantile(y, 0.95)) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f0(tt)+bias.0, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')

#confidence from paraboot:
sapply(coefs.boot, function(x) x[,2]) %>% range -> yy
sapply(coefs.boot, function(x) x[,2]) %>% apply(1, function(y) quantile(y, 0.05)) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[1]), ylim=yy, bty='n')
par(new=TRUE)
sapply(coefs.boot, function(x) x[,2]) %>% apply(1, function(y) quantile(y, 0.95)) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f1(tt)+bias.1, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')

#confidence from paraboot:
sapply(coefs.boot, function(x) x[,3]) %>% range -> yy
sapply(coefs.boot, function(x) x[,3]) %>% apply(1, function(y) quantile(y, 0.05)) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[2]), ylim=yy, bty='n')
par(new=TRUE)
sapply(coefs.boot, function(x) x[,3]) %>% apply(1, function(y) quantile(y, 0.95)) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f2(tt)+bias.2, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')

#dev.off()



## @knitr linked-bootstrap-plot-efron

#confidence from Efron:
sapply(coefs.boot, function(x) x[,1]) %>% range -> yy
(sapply(m$fits, function(x) x$coef)[1,] - 1.96*sapply(sd.boot, function(x) x[1,1])) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[0]), ylim=yy, bty='n')
par(new=TRUE)
(sapply(m$fits, function(x) x$coef)[1,] + 1.96*sapply(sd.boot, function(x) x[1,1])) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f0(tt)+bias.0, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')


#confidence from Efron:
sapply(coefs.boot, function(x) x[,2]) %>% range -> yy
(sapply(m$fits, function(x) x$coef)[2,] - 1.96*sapply(sd.boot, function(x) x[2,2])) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[1]), ylim=yy, bty='n')
par(new=TRUE)
(sapply(m$fits, function(x) x$coef)[2,] + 1.96*sapply(sd.boot, function(x) x[2,2])) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f1(tt)+bias.1, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')


#confidence from Efron:
sapply(coefs.boot, function(x) x[,3]) %>% range -> yy
(sapply(m$fits, function(x) x$coef)[3,] - 1.96*sapply(sd.boot, function(x) x[3,3])) %>% plot(x=tt, type='l', xlab='t', ylab=expression(beta[2]), ylim=yy, bty='n')
par(new=TRUE)
(sapply(coefs.boot, function(x) x[,3]) %>% rowMeans + 1.96*sapply(sd.boot, function(x) x[3,3])) %>% plot(x=tt, type='l', ylim=yy, ann=FALSE, bty='n', xaxt='n', yaxt='n')
par(new=TRUE)
plot(x=tt, y=f2(tt)+bias.2, ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n')





## @knitr plot-bias
layout(matrix(1:3,1,3))

#intercept plot:
yy = range(f0(tt))
plot(x=tt, y=f0(tt), type='l', xlab='t', ylab='Intercept', bty='n', ylim=yy)
par(new=TRUE)
plot(x=tt, y=f0(tt) + bias.0, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='blue', lwd=2, lty=2)
par(new=TRUE)
plot(x=tt, y=f0(tt) + bias.0 + bias.0.2, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='red', lwd=2, lty=3)
legend(x='topleft', c("Truth", "First order", "Second order"), lty=c(1,2,3), lwd=c(1,2,2), col=c('black', 'blue', 'red'), bty='n')

#beta_1 plot:
yy = range(f1(tt))
plot(x=tt, y=f1(tt), type='l', xlab='t', ylab='Intercept', bty='n', ylim=yy)
par(new=TRUE)
plot(x=tt, y=f1(tt) + bias.1, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='blue', lwd=2, lty=2)
par(new=TRUE)
plot(x=tt, y=f1(tt) + bias.1 + bias.1.2, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='red', lwd=2, lty=3)
legend(x='topleft', c("Truth", "First order", "Second order"), lty=c(1,2,3), lwd=c(1,2,2), col=c('black', 'blue', 'red'), bty='n')

#beta_2 plot:
yy = range(f2(tt))
plot(x=tt, y=f2(tt), type='l', xlab='t', ylab='Intercept', bty='n', ylim=yy)
par(new=TRUE)
plot(x=tt, y=f2(tt) + bias.2, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='blue', lwd=2, lty=2)
par(new=TRUE)
plot(x=tt, y=f2(tt) + bias.2 + bias.2.2, type='l', ann=FALSE, xaxt='n', yaxt='n', bty='n', ylim=yy, col='red', lwd=2, lty=3)
legend(x='topright', c("Truth", "First order", "Second order"), lty=c(1,2,3), lwd=c(1,2,2), col=c('black', 'blue', 'red'), bty='n')

# 
# ## @knitr scratch
# 
# #Plot draws from the linked bootstrap: TOO SMOOTH
# yy2 = sapply(beta.star.1, function(x) range(x[2,])) %>% range
# plot(x=tt,y=f1(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy2, type='l')
# par(new=TRUE)
# plot(x=tt, y=coefs[[B]][,2], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue')
# par(new=TRUE)
# plot(x=tt, y=f1(tt) + bias.1, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue', lty=2)
# for (b in 1:B) {
#     par(new=TRUE)
#     plot(x=tt, y=beta.star.1[[b]][2,], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='red')
# }
# 
# 
# 
# #Now a full distribution bootstrap (too rough?):
# Sigma = list()
# for (i in 1:n) {
#     Sigma[[i]] = with(summary(m[['fits']][[i]][['model']][['adamodel']]), cov.unscaled * dispersion)[1:3,1:3]
# }
# 
# SS = bdiag(Sigma)
# for (i in 1:(n-1)) {
#     X.i = cbind(1, X1, X2)
#     
#     for (j in (i+1):n) {
#         X.j = cbind(1, X1, X2)
#         scale = sqrt(summary(m[['fits']][[i]][['model']][['adamodel']])$dispersion *
#                          summary(m[['fits']][[j]][['model']][['adamodel']])$dispersion)
#         
#         unscaled = summary(m[['fits']][[i]][['model']][['adamodel']])$cov.unscaled[1:3,1:3] %*%
#             t(X.i) %*% sqrt(W[[i]] %*% W[[j]]) %*% 
#             X.j %*% summary(m[['fits']][[j]][['model']][['adamodel']])$cov.unscaled[1:3,1:3]
#         
#         SS[((i-1)*3+1):(3*i), ((j-1)*3+1):(3*j)] =
#             SS[((j-1)*3+1):(3*j), ((i-1)*3+1):(3*i)] = scale * unscaled
#     }
# }
# 
# SS.posdef = SS %*% t(SS)
# SS.eigen = eigen(SS.posdef)
# 
# sqSS = SS.eigen$vectors %*% diag(SS.eigen$values^(1/4)) %*% t(SS.eigen$vectors)
# 
# #Draw from a joint parametric distribution:
# beta.star.3 = list()
# beta.star.4 = list()
# for (k in 1:B) {
#     #beta.star.3[[k]] = (t(coefs) + matrix(sqSS %*% rep(rnorm(6),n), 6, n))
#     beta.star.4[[k]] = (t(coefs[[B]][,1:3]) + matrix(sqSS %*% rnorm(3*n), 3, n))
# }
# 
# yy2 = sapply(beta.star.4, function(x) range(x[1,])) %>% range
# plot(x=tt,y=f0(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy2, type='l')
# par(new=TRUE)
# plot(x=tt, y=coefs[[B]][,1], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='blue')
# for (b in 1:B) {
#     par(new=TRUE)
#     plot(x=tt, y=beta.star.4[[b]][1,], bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy2, type='l', col='red')
# }
# 
# 
# 
# g1 = gam(resid[[B]]~s(tt, k=50))
# s2 = g1$sig2
# 
# plot(x=tt, y=resid[[B]], bty='n', xlab='t', ylab='residual')
# par(new=TRUE)
# plot(x=tt, y=fitted(g1), type='l', ylim=range(resid), ann=FALSE, xaxt='n', yaxt='n', bty='n')
