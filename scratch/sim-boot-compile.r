library(RandomFields)
library(dplyr)
library(lagr)
library(doMC)

registerDoMC(6)

N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)

setwd("~/Desktop/sim-out/")
files = dir()
runs = (1:12)[-2]

for (s in runs) {
    set.seed(s)
    
    #Calculate the coefficient surfaces
    B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
    B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
    B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]
    
    #s.files = grep(paste("-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""), files, perl=TRUE)
    res.files = grep(paste("residual-bootstrap-coefs\\.", s, "\\.[0-9]+\\.csv", sep=""),
                     files, perl=TRUE)
    
    ff = gregexpr(paste("residual-bootstrap-coefs\\.", s, "\\.(?<bootstrap>[0-9]+)\\.csv", sep=""),
                     files[res.files], perl=TRUE)
    
    bootstraps = sapply(1:length(res.files), function(k) 
        substr(files[res.files[k]], attr(ff[[k]], 'capture.start')[1,'bootstrap'],
               attr(ff[[k]], 'capture.start')[1,'bootstrap'] +
                   attr(ff[[k]], 'capture.length')[1,'bootstrap'] - 1)) %>% as.numeric
    
    res.boot = array(NA, dim=c(400, 5, length(res.files)))
    for (b in 1:length(res.files)) {
        res.boot[,,b] = as.matrix(read.table(files[res.files[b]], header=TRUE, row.names=1))
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
    
    bayes.boot = array(NA, dim=c(400, 5, length(bayes.files)))
    for (b in 1:length(bayes.files)) {
        bayes.boot[,,b] = as.matrix(read.table(files[bayes.files[b]], header=TRUE, row.names=1))
    }
    
    
    
    
    
    
    #BIAS CORRECTION
    load(paste("bootstrap-junk.", s, ".Rdata", sep=""))
    bw = trace$bw[which.min(trace$loss)]
    
    raw.data = read.table(paste("raw-data.", s, ".csv", sep=""), header=TRUE, row.names=1)
    raw.coef = read.table(paste("raw-model-coefs.", s, ".csv", sep=""), header=TRUE, row.names=1)
    
    raw.model = lagr(Y~X1+X2+X3+X4, 
                                 data=raw.data, 
                                 family='gaussian', 
                                 coords=c('loc.x','loc.y'), 
                                 longlat=FALSE,
                                 varselect.method='AICc', 
                                 bw=bw,
                                 bw.type='knn',
                                 kernel=epanechnikov,
                                 verbose=FALSE, 
                                 n.lambda=100, 
                                 lagr.convergence.tol=0.005)
    
    
    #Empirical second derivative of beta.0
    slx = sapply(raw.model$fits, function(x) x$model$adamodel$coef[8])
    sly = sapply(raw.model$fits, function(x) x$model$adamodel$coef[9])
    cc.1 = vector()
    for (i in 1:N**2) {
        h = raw.model$fits[[i]]$bw
        X = cbind(raw.data$loc.x - raw.data$loc.x[i], raw.data$loc.y - raw.data$loc.y[i])
        w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
        
        m2x = lm(slx~X, weights=w)
        m2y = lm(sly~X, weights=w)
        cc.1 = c(cc.1, m2x$coef[2] + m2y$coef[3])
    }
    
    empirical.bias.1 = cc.1 * sapply(raw.model$fits, function(x) x$bw**2) / 10
    
    
    #Confidence intervals:
    res.CI2 = sapply(1:400, function(k) quantile(res.boot[k,2,], probs=c(0.05, 0.95)))
    bayes.CI2 = sapply(1:400, function(k) quantile(bayes.boot[k,2,], probs=c(0.05, 0.95)))
    
    
    
    #Confidence intervals:
    res.CI = sapply(1:400, function(k) quantile(res.boot[k,2,], probs=c(0.05, 0.95)) - empirical.bias.1[k])
    bayes.CI = sapply(1:400, function(k) quantile(bayes.boot[k,2,], probs=c(0.05, 0.95)) - empirical.bias.1[k])
}