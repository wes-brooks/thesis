library(expm)
library(lagr)
library(doMC)
library(MASS)

registerDoMC(7)

source("scratch/inference/simulate-data.r")

bw = lagr.tune(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bwselect.method="AIC", kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005, tol.bw=0.0005)
bw$trace[order(bw$trace[,1]),][,c(1,2)] -> trace
bw.range = range(trace[trace[,2]<min(trace[,2])+10,1])

#Number of bootstrap draws:
B = 25

bootstrap.data = sim
hh = seq(bw.range[1], bw.range[2], length.out=10)
for (h in hh) {
    cc = list()
    model = lagr(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005)
    
    bootstrap.seed = sapply(1:B, function(m) return(matrix(rnorm(15, 0, 1))))
    bootstrap.response = matrix(NA, 0, B)
    for (i in 1:nrow(sim)) {
        mu = model[['model']][[i]][['model']][['LS.coefs']][1:5]
        Sigma = 4/3/pi * model[['model']][[i]][['model']][['results']][['full.model.cov']]
        S = sqrtm(Sigma)
        bootstrap.coefs = t((S[1:5,] %*% bootstrap.seed) + mu)
        cc[[i]] = bootstrap.coefs
        eta = bootstrap.coefs %*% matrix(as.numeric(cbind(1, sim[i,2:5])))
        bootstrap.response = rbind(bootstrap.response,
            rnorm(n=length(eta), mean=eta, sd=sqrt(model[['model']][[i]][['dispersion']])))
    }

    bootstrap.model = list()
    fl = sim[133,c('loc.x', 'loc.y')]
    for (i in 1:B) {
        cat(paste(i, "\n", sep=""))
        bootstrap.data$Y = bootstrap.response[,i]
        bootstrap.model[[i]] = lagr(Y~X1+X2+X3+X4, data=bootstrap.data, family='gaussian', coords=c('loc.x','loc.y'), fit.loc=fl, longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005)
    }
    
    is.zero = matrix(NA,5,0)
    coef.distribution = matrix(NA,5,0)
    for (j in 1:B) {
        
        coefs = bootstrap.model[[j]][['model']][[1]][['coef']]
        bootstrap.model[[j]][['model']][[1]][['tunelist']][['criterion']] -> crit
        is.zero=cbind(is.zero, matrix(as.double(coefs==0), nrow(coefs), ncol(coefs)) %*% exp(0.5*(min(crit)-crit)) / sum(exp(0.5*(min(crit)-crit))))
        coef.distribution = cbind(coef.distribution, as.double(coefs %*% exp(0.5*(min(crit)-crit)) / sum(exp(0.5*(min(crit)-crit)))))
    }
    
    #Following Efron (2014):
    V = 4/3/pi * model[['model']][[133]][['model']][['results']][['full.model.cov']][1:5,1:5]
    cov(t(coef.distribution), cc[[133]])
}
