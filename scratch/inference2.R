library(lagr)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(8)

#Generate the covariates:
N = 14 # number of width and length divisions in the domain
B = 200 # number of bootstrap resamplings of the original data
S = 80 # steps in the grid search for bandwidth
coord = seq(0, 1, length.out=N)

# #Seed the RNG
set.seed(11181982)

#Generate the bootstrap draws:
index = list()
index[[1]] = 1:N**2
if (B > 1) {
    for (i in 2:B) {
        index[[i]] = sample(1:N**2, replace=TRUE)
    }
}

if (interactive()) {
    i=1
    j=1
} else {
    iter = commandArgs(trailingOnly=TRUE)[1] 
    print(iter)
    iter = as.numeric(iter)
    i = ((iter-1) %/% S) + 1
    j = ((iter-1) %% S) + 1
    
    print(i)
    print(j)
}

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
grid = cbind(loc.x, loc.y)

sim = read.table("~/git/gwr/scratch/sim.txt")

hh = c(0.1726461, 0.1830205, 0.2284030, 0.2084279, 0.1985406, 0.1871199, 0.2153826, 0.1800163, 0.2327630, 0.2081909)


models = list()

for (k in 1:10) {
    h = hh[k]
    print(h)
    indx = index[[i]]
    data = sim[indx,]
    models[[i]] = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), fit.loc=grid, longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005, jacknife=FALSE)

#df = sum(sapply(models[[i]][['model']], function(x) tail(x[['tunelist']][['df-local']], 1)))
df2 = sum(sapply(1:length(models[[i]][['model']]), function(k) models[[i]][['model']][[k]][['tunelist']][['df']]*sum(indx==k)/models[[i]][['model']][[k]][['weightsum']]))

#Write LAGR coefficients:
coefs = t(sapply(models[[i]][['model']], function(x) x[['coef']]))

#Write LAGR results:
#fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, data[k,2:5])))
fitted2 = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, sim[k,2:5])))
fitted3 = fitted2[unique(indx)]

#dev.resids = gaussian()$dev.resids(data$Y, fitted, rep(1,N**2))
dev.resids2 = gaussian()$dev.resids(sim$Y, fitted2, rep(1,N**2))
dev.resids3 = gaussian()$dev.resids(sim$Y[unique(indx)], fitted3, w)

#ll = gaussian()$aic(data$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))
ll2 = gaussian()$aic(sim$Y, N**2, fitted2, rep(1,N**2), sum(dev.resids2))
ll3 = gaussian()$aic(sim$Y[unique(indx)], length(unique(indx)), fitted3, w, sum(dev.resids3))

#err = ll
aic = c(aic, ll2 + 2*df2)
print(ll)

# write(c(h, err), "~/git/gwr/scratch/trace.txt", append=TRUE)


}