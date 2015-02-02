library(RandomFields)
library(dplyr)
library(lagr)
library(doMC)

registerDoMC(3)

N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)
x.c = rep(coord, each=N)
y.c = rep(coord, times=N)
location = data.frame(x=x.c, y=y.c)
n = N**2

setwd("~/Desktop/sim-out/")
files = dir()
runs = 1:5

res.boot = list()
bayes.boot = list()
orc.boot = list()

#Calculate the coefficient surfaces
B1.fun = function(x,y) 3.75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
    3.75*exp(-(9*x+1)^2/49 - (9*y+1)/10) +
    2.5*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
    exp(-(9*x-4)^2 - (9*y-7)^2) - 0.5

B2.fun = function(x,y) 3*sort(c(0, 5*((x-3/4)^2+(y-1/4)^2) -1/2,-5*((x-1/4)^2 + (y-3/4)^2) + 1/2 ))[2]

B1.. = function(x,y) -3.75 * 18/4 * (9 - 18/4*(9*x-2)^2) * exp(-((9*x-2)^2 + (9*y-2)^2)/4) -
    3.75 * 18/4 * (9 - 18/4*(9*y-2)^2) * exp(-((9*x-2)^2 + (9*y-2)^2)/4) -
    
    3.75 * 18/49 * (9 - 18/49*(9*x+1)^2) * exp(-(9*x+1)^2/49 - (9*y+1)/10) +
    3.75 * 81/100 * exp(-(9*+1)^2/49 - (9*y+1)/10) -
    
    2.5 * 18/4 * (9 - 18/4*(9*x-7)^2) * exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
    2.5 * 18/4 * (9 - 18/4*(9*y-3)^2) * exp(-((9*x-7)^2 + (9*y-3)^2)/4) +
    
    18 * (9 - 18*(9*x-4)^2) * exp(-(9*x-4)^2 - (9*y-7)^2) +
    18 * (9 - 18*(9*y-7)^2) * exp(-(9*x-4)^2 - (9*y-7)^2)

B1 = B1.fun(location$x, location$y)
B2 = sapply(1:nrow(location), function(k) B2.fun(location$x[k], location$y[k]))
B3 = rep(-0.1, n)

bias.0 = B1..(location$x, location$y) * sapply(model[['fits']], function(x) (x[['bw']]*2)**2) / 10

for (s in runs) {
    set.seed(s)
    
    res.files = grep(paste("residual-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""),
                     files, perl=TRUE)
    
    ff = gregexpr(paste("residual-bootstrap-coefs\\.", s, "\\.(?<bootstrap>[0-9]+)\\.csv", sep=""),
                     files[res.files], perl=TRUE)
    
    bootstraps = sapply(1:length(res.files), function(k) 
        substr(files[res.files[k]], attr(ff[[k]], 'capture.start')[1,'bootstrap'],
               attr(ff[[k]], 'capture.start')[1,'bootstrap'] +
                   attr(ff[[k]], 'capture.length')[1,'bootstrap'] - 1)) %>% as.numeric
    
    res.boot[[s]] = array(NA, dim=c(400, 5, length(res.files)))
    for (b in 1:length(res.files)) {
        res.boot[[s]][,,b] = tryCatch(
            as.matrix(read.table(files[res.files[b]], header=TRUE, row.names=1)),
            error = function(e) matrix(NA, 400, 5)
        )
    }
    
    
    
    
    
    #s.files = grep(paste("-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""), files, perl=TRUE)
    bayes.files = grep(paste("bayesian-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""),
                     files, perl=TRUE)
    
    ff = gregexpr(paste("bayesian-bootstrap-coefs\\.", s, "\\.(?<bootstrap>[0-9]+)\\.csv", sep=""),
                  files[bayes.files], perl=TRUE)
    
    bootstraps = sapply(1:length(bayes.files), function(k) 
        substr(files[bayes.files[k]], attr(ff[[k]], 'capture.start')[1,'bootstrap'],
               attr(ff[[k]], 'capture.start')[1,'bootstrap'] +
                   attr(ff[[k]], 'capture.length')[1,'bootstrap'] - 1)) %>% as.numeric
    
    bayes.boot[[s]] = array(NA, dim=c(400, 5, length(bayes.files)))
    for (b in 1:length(bayes.files)) {
        bayes.boot[[s]][,,b] = tryCatch(
            as.matrix(read.table(files[bayes.files[b]], header=TRUE, row.names=1)),
            error = function(e) matrix(NA, 400, 5)
        )
    }
    
    
    
    
    
    #s.files = grep(paste("-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""), files, perl=TRUE)
    orc.files = grep(paste("oracle-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""),
                     files, perl=TRUE)
    
    ff = gregexpr(paste("oracle-bootstrap-coefs\\.", s, "\\.(?<bootstrap>[0-9]+)\\.csv", sep=""),
                  files[orc.files], perl=TRUE)
    
    bootstraps = sapply(1:length(orc.files), function(k) 
        substr(files[orc.files[k]], attr(ff[[k]], 'capture.start')[1,'bootstrap'],
               attr(ff[[k]], 'capture.start')[1,'bootstrap'] +
                   attr(ff[[k]], 'capture.length')[1,'bootstrap'] - 1)) %>% as.numeric
    
    orc.boot[[s]] = array(NA, dim=c(400, 5, length(res.files)))
    for (b in 1:length(res.files)) {
        orc.boot[[s]][,,b] = tryCatch(
            as.matrix(read.table(files[orc.files[b]], header=TRUE, row.names=1)),
            error = function(e) matrix(NA, 400, 5)
        )
    }
    
    
    
    
    
#     #BIAS CORRECTION
#     load(paste("bootstrap-junk.", s, ".Rdata", sep=""))
#     bw = trace$bw[which.min(trace$loss)]
#     
#     raw.data = read.table(paste("raw-data.", s, ".csv", sep=""), header=TRUE, row.names=1)
#     raw.coef = read.table(paste("raw-model-coefs.", s, ".csv", sep=""), header=TRUE, row.names=1)
#     
#     raw.model = lagr(Y~X1+X2+X3+X4, 
#                                  data=raw.data, 
#                                  family='gaussian', 
#                                  coords=c('loc.x','loc.y'), 
#                                  longlat=FALSE,
#                                  varselect.method='AICc', 
#                                  bw=bw,
#                                  bw.type='knn',
#                                  kernel=epanechnikov,
#                                  verbose=FALSE, 
#                                  n.lambda=100, 
#                                  lagr.convergence.tol=0.005)
#     
#     
cc = model$fits %>% sapply(function(x) x$coef) %>% t

    #Empirical second derivative of beta.0
    slx = sapply(model$fits, function(x) x$model$adamodel$coef[8])
    sly = sapply(model$fits, function(x) x$model$adamodel$coef[9])
    cc.0 = vector()
    for (i in 1:N**2) {
        h = model$fits[[i]]$bw
        X = cbind(location$x - location$x[i], location$y - location$y[i])
        w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
        
        m2x = lm(slx~X, weights=w)
        m2y = lm(sly~X, weights=w)
        cc.0 = c(cc.0, m2x$coef[2] + m2y$coef[3])
    }
    empirical.bias.0 = cc.0 * sapply(model$fits, function(x) x$bw**2) / 10

cc.1 = vector()
for (i in 1:N**2) {
    h = model$fits[[i]]$bw / 2
    X = cbind(location$x - location$x[i], location$y - location$y[i])
    w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
    
    m2x = lm(slx~X, weights=w)
    m2y = lm(sly~X, weights=w)
    cc.1 = c(cc.1, m2x$coef[2] + m2y$coef[3])
}
empirical.bias.1 = cc.1 * sapply(model$fits, function(x) x$bw**2) / 10


cc.2 = vector()
for (i in 1:N**2) {
    h = model$fits[[i]]$bw / 4
    X = cbind(location$x - location$x[i], location$y - location$y[i])
    w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
    
    m2x = lm(slx~X, weights=w)
    m2y = lm(sly~X, weights=w)
    cc.2 = c(cc.2, m2x$coef[2] + m2y$coef[3])
}
empirical.bias.2 = cc.2 * sapply(model$fits, function(x) x$bw**2) / 10

pdf(paste("~/Desktop/bias.", s, ".pdf", sep=""), width=18,height=6)
empirical.bias.0 %>% plot(type='l', col='red', ylim=c(-2,2))
par(new=TRUE)
empirical.bias.1 %>% plot(type='l', col='blue', ylim=c(-2,2))
par(new=TRUE)
empirical.bias.2 %>% plot(type='l', col='green', ylim=c(-2,2))
par(new=TRUE)
(cc[,2]-B1) %>% plot(type='l', ylim=c(-2,2))
dev.off()
    
    #Confidence intervals:
    res.CI2 = sapply(1:400, function(k) quantile(res.boot[k,2,], probs=c(0.05, 0.95)))
    bayes.CI2 = sapply(1:400, function(k) quantile(bayes.boot[k,2,], probs=c(0.05, 0.95)))
    
    
    
    #Confidence intervals:
    res.CI = sapply(1:400, function(k) quantile(res.boot[k,2,], probs=c(0.05, 0.95)) - empirical.bias.1[k])
    bayes.CI = sapply(1:400, function(k) quantile(bayes.boot[k,2,], probs=c(0.05, 0.95)) - empirical.bias.1[k])
}




