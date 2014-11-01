library(expm)

bw = lagr.tune(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bwselect.method="AIC", kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005, tol.bw=0.0005)
bw$trace[order(bw$trace[,1]),][,c(1,2)] -> trace
bw.range = range(trace[trace[,2]<min(trace[,2])+10,1])


bootstrap.data = sim
hh = seq(bw.range[1], bw.range[2], length.out=10)
for (h in hh) {
    
    model = lagr(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005)
    bootstrap.response = t(sapply(1:nrow(sim), function(i) {
        coefs = mvrnorm(20, mu=model[['model']][[i]][['model']][['LS.coefs']][1:5], Sigma=4/3/pi * model[['model']][[i]][['model']][['results']][['full.model.cov']][1:5,1:5])
        eta = coefs %*% matrix(as.numeric(cbind(1, sim[i,2:5])))
        rnorm(n=length(eta), mean=eta, sd=sqrt(model[['model']][[i]][['dispersion']]))
    }))

    bootstrap.model = list()
    for (i in 1:20) {
        cat(paste(i, "\n", sep=""))
        bootstrap.data$Y = bootstrap.response[,i]
        bootstrap.model[[i]] = lagr(Y~X1+X2+X3+X4, data=bootstrap.data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005)
    }
    
    is.zero = matrix(NA,5,0)
    coef.distribution = matrix(NA,5,0)
    for (j in 1:20) {
        
        coefs = bootstrap.model[[j]][['model']][[1]][['coef']]
        bootstrap.model[[j]][['model']][[1]][['tunelist']][['criterion']] -> crit
        is.zero=cbind(is.zero, matrix(as.double(coefs==0), nrow(coefs), ncol(coefs)) %*% exp(0.5*(min(crit)-crit)) / sum(exp(0.5*(min(crit)-crit))))
        coef.distribution = cbind(coef.distribution, as.double(coefs %*% exp(0.5*(min(crit)-crit)) / sum(exp(0.5*(min(crit)-crit)))))
    }
}
