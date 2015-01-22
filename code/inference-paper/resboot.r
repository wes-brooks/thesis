## @knitr simulate-data
#Set constants:
B1 = 1
n = 100
tt = seq(0, 5, len=n)

#Coefficient functions:
f0 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f1 = function(x) {0.2*x^2 - x + cos(x)}
f2 = function(x) {(x+1)^(-1) - 1/6}

f0.second.derivative = function(x) {-2.4*x^2 + 6*x + 1.4}
f1.second.derivative = function(x) {0.4 - cos(x)}
f2.second.derivative = function(x) {2*(x+1)^(-3)}

#Generate the raw data sets (not bootstrap)
set.seed(11181982)

#simulate data
X1 = rnorm(n, mean=3, sd=2)
X2 = rnorm(n, mean=3, sd=2)
y = f0(tt) + X1*f1(tt) + X2*f2(tt) + rnorm(n, mean=0, sd=1+tt/10)
df = data.frame(y, x1=X1, x2=X2, t=tt)

bw = lagr.tune(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc',
                    kernel=epanechnikov, bw.type='dist', bwselect.method='AIC', tol.bw=0.1, verbose=FALSE,
                    lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)

#Fit a VCR model to the simulated data by LAGR
m = lagr(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc',
         bw=bw, kernel=epanechnikov, bw.type='dist', verbose=TRUE,
         lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)

#Get the observation weights for each local fit in the VCR model:
W = list()
for (i in 1:n) {
    h = m[['fits']][[i]][['bw']]
    w = epanechnikov(abs(tt-tt[i]), h)
    W[[i]] = diag(w)
}

#Store the fitted values, residuals, and coefficient estimates:
fitted = sapply(m[['fits']], function(x) tail(x[['fitted']],1))
resid = df$y - fitted
coefs = t(sapply(m[['fits']], function(x) x[['model']][['beta']][,ncol(x[['model']][['beta']])]))

## @knitr asymptotic-bias

#Compute asymptotic bias from the true coefficient functions:
bias.0 = f0.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.1 = f1.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10
bias.2 = f2.second.derivative(tt) * sapply(m[['fits']], function(x) x[['bw']]**2) / 10







## @knitr empirical-bias

#Empirical bias:
empirical.bias.0 = list()
empirical.bias.1 = list()
empirical.bias.2 = list()
empirical.bias = list()

#Empirical second derivative of beta.0
sl = coefs[,4]
cc.0 = matrix(0,0,2)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = tt - tt[i]
    
    m2 = lm(sl~X, weights=w)
    cc.0 = rbind(cc.0, m2$coef)
}

empirical.bias.0 = cc.0[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10

#Empirical second derivative of beta.1
sl = coefs[,5]
cc.1 = matrix(0,0,2)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = tt - tt[i]
    
    m2 = lm(sl~X, weights=w)
    cc.1 = rbind(cc.1, m2$coef)
}

empirical.bias.1 = cc.1[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10

#Empirical second derivative of beta.2
sl = coefs[,6]
cc.2 = matrix(0,0,2)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = tt - tt[i]
    
    m2 = lm(sl~X, weights=w)
    cc.2 = rbind(cc.2, m2$coef)
}

empirical.bias.2 = cc.2[,2] * sapply(m[['fits']], function(x) x[['bw']]**2) / 10

empirical.bias = cbind(empirical.bias.0, empirical.bias.1, empirical.bias.2)





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
legend(x='bottomleft', legend=c("Truth", "Estimable", "Estimate", "Bias corrected"), lty=c(1,2,1,2), col=c("black", "blue", "red", "green"), lwd=c(1,2,1,2), bty='n')

#Plot a realization of the fitted model: \beta_1
yy = range(c(f1(tt), sapply(coefs, function(x) x[,2])))
plot(x=tt, y=f1(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,2])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
par(new=TRUE)
plot(x=tt, y=f1(tt) + bias.1, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,2])) - rowMeans(sapply(empirical.bias, function(b) b[,2])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='green', lty=2, lwd=2)
legend(x='topright', legend=c("Truth", "Estimable", "Estimate", "Bias corrected"), lty=c(1,2,1,2), col=c("black", "blue", "red", "green"), lwd=c(1,2,1,2), bty='n')

#Plot a realization of the fitted model: \beta_2
yy = range(c(f2(tt), sapply(coefs, function(x) x[,3])))
plot(x=tt,y=f2(tt), bty='n', xlab='t', ylab='Intercept', ylim=yy, type='l')
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,3])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='red')
par(new=TRUE)
plot(x=tt, y=f2(tt) + bias.2, bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='blue', lty=2, lwd=2)
par(new=TRUE)
plot(x=tt, y=rowMeans(sapply(coefs, function(c) c[,3])) - rowMeans(sapply(empirical.bias, function(b) b[,3])), bty='n', ann=FALSE, yaxt='n', xaxt='n', ylim=yy, type='l', col='green', lty=2, lwd=2)
legend(x='topright', legend=c("Truth", "Estimable", "Estimate", "Bias corrected"), lty=c(1,2,1,2), col=c("black", "blue", "red", "green"), lwd=c(1,2,1,2), bty='n')





## @knitr residual-bootstrap-sample
B = 200

#Linked bootstrap draws to generate the resampled response:
Y.star = matrix(0, 0, B)

for (j in 1:n) {
    cat(paste("j is ", j, "\n", sep=""))
    Y.star = rbind(Y.star, m$fits[[j]]$model$adamodel$fitted[as.character(j)] + 
                       sample(m$fits[[j]]$model$adamodel$residuals, B, replace=TRUE))
}


AIC.boot = vector()
AICc.boot = vector()
df.boot = vector()

#Compute a bandwidth distribution:
b=1
indx = which(bw$trace$loss < (min(bw$trace$loss+20)))
max.lik = as.data.frame(bw$trace[indx,1:2])
m.bw = lm(loss ~ log(bw) + I(log(bw)^2), data=max.lik)
mu = -m.bw$coef[2] / m.bw$coef[3] / 2
sd = sqrt(1 / m.bw$coef[3])
bw.boot = bw.b = exp(rnorm(B, mean=mu, sd=sd))

#Run estimation on the parametric bootstrap draws:
coefs.boot = list()
conf.zero = list()
for (b in 30:90) {   
    print(b) 
    print(bw.boot[b])
    df.b = df
    df.b$y = Y.star[,b]
    m.b = lagr(y~x1+x2, data=df.b, family='gaussian', coords='t', varselect.method='wAICc', kernel=epanechnikov, bw.type='dist', bw=bw.boot[b], verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    
    AIC.boot = c(AIC.boot, m.b$AIC)
    AICc.boot = c(AICc.boot, m.b$AICc)
    df.boot = c(df.boot, m.b$df)
    
    coefs.boot[[b]] = t(sapply(m.b$fits, function(x) x$coef))
    conf.zero[[b]] = t(sapply(m.b$fits, function(x) x$conf.zero))
}




#sd by Efron method (nonparametric delta method) (nonparametric version):
sd.boot = list()
sd.zero.boot = list()

for (k in 1:n) {
    t.b = sapply(coefs.boot, function(x) x[k,]) %>% t
    t.mean = matrix(rep(colMeans(t.b),each=B),B,3)
    
    z.b = sapply(conf.zero, function(x) x[k,]) %>% t
    z.mean = matrix(rep(colMeans(z.b),each=B),B,3)
    
    cv = matrix(0,3,3)
    cz = matrix(0,3,3)
    
    for (i in 1:n) {        
        Y.b = sapply(wt, function(x) x[i])
        Y.b = Y.b - mean(Y.b)
        
        cov.b = t(Y.b) %*% (t.b - t.mean) / B
        coz.b = t(Y.b) %*% (z.b - z.mean) / B
        
        cv = cv + t(cov.b) %*% cov.b
        cz = cz + t(coz.b) %*% coz.b
    }
    
    sd.boot[[k]] = sqrtm(cv)
    sd.zero.boot[[k]] = sqrtm(cz)  
}



## @knitr bw-boot-Dirichlet
#Run estimation on the parametric bootstrap draws:
interval = 10
bw.trace.boot = list()
bw.range = c(0,2)
for (b in 1:(B %/% interval)) {   
    print(b) 
    
    w = wt[[interval*b]]
    
    bw.trace.boot[[b]] = lagr.tune(y~x1+x2, data=df, weights=w,
                                   family='gaussian', coords='t',
                                   varselect.method='wAICc',
                                   kernel=epanechnikov, bw.type='dist',
                                   bwselect.method='AIC', verbose=TRUE,
                                   lagr.convergence.tol=0.005,
                                   lambda.min.ratio=0.01,
                                   n.lambda=80, range=bw.range)
    
    tr = bw.trace.boot[[b]]$trace
    bw.range = c(0, max(2, max(tr$bw[tr$loss < min(tr$loss)+20])))
}

sapply(bw.boot, function(x) sum(dnorm(bw.trace.boot %>% sapply(function(y) y$bw), mean=x, sd=kdest$bw))/20) / (bw.boot %>% sapply(function(x) dnorm(log(x) , mean=mu, sd=sd))) -> importance




## @knitr bayesian-bootstrap-plot-efron
layout(matrix(1:3, 1, 3))
legend.loc = c('bottomleft', 'topright', 'topright')
f = c(f0, f1, f2)
b = list(bias.0, bias.1, bias.2)
ylabs = c(expression(beta[0]), expression(beta[1]), expression(beta[2]))
#confidence from Efron:
for (k in 1:3) {
    sapply(coefs.boot, function(x) x[,k]) %>% range -> yy
    #plot(rowMeans(sapply(coefs.boot, function(x) x[,1])), x=tt, type='l', xlab='t', ylab=expression(beta[0]), ylim=yy, bty='n')
    #par(new=TRUE)
    (rowMeans(sapply(coefs.boot, function(x) x[,k])) - 1.96*sapply(sd.boot, function(x) diag(x)[k])) %>% plot(x=tt, type='l', ylim=yy, bty='n', xlab='t', ylab=ylabs[k], lty=1)
    par(new=TRUE)
    (rowMeans(sapply(coefs.boot, function(x) x[,k])) + 1.96*sapply(sd.boot, function(x) diag(x)[k])) %>% plot(x=tt, type='l', ylim=yy, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=1)
    par(new=TRUE)
    plot(x=tt, y=f[[k]](tt)+b[[k]], ylim=yy, type='l', col='red', ann=FALSE, bty='n', xaxt='n', yaxt='n', lty=3, lwd=2)
    par(new=TRUE)
    plot(x=tt, y=f[[k]](tt), ylim=yy, type='l', col='blue', ann=FALSE, bty='n', xaxt='n', yaxt='n', lty=2, lwd=2)
    legend(c("95% confidence band", "truth", "biased"), x=legend.loc[k], col=c("black", "blue", "red"), lty=c(1,2,3), lwd=c(1,2,2), bty='n')
}



## @knitr scratch

#sd by Efron method (nonparametric delta method) (parametric version):
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





## @knitr linked-bootstrap-estimate

#Run estimation on the parametric bootstrap draws:
coefs.boot.par = list()
conf.zero.par = list()
df.b = df
for (b in 1:B) {
    print(b)
    df.b$y = Y.star.1[[b]]
    
    m.b = lagr(y~x1+x2, data=df.b, family='gaussian', coords='t', varselect.method='wAICc', kernel=epanechnikov, bw.type='dist', bw=bw.boot[b], verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)
    coefs.boot.par[[b]] = t(sapply(m.b$fits, function(x) x$coef))
    conf.zero.par[[b]] = t(sapply(m.b$fits, function(x) x$conf.zero))
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

