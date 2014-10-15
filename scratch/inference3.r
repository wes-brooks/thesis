library(lagr)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(8)

#Generate the covariates:
N = 14 # number of width and length divisions in the domain
B = 20 # number of bootstrap resamplings of the original data
S = 31 # steps in the grid search for bandwidth
coord = seq(0, 1, length.out=N)

# #Seed the RNG
set.seed(111882)

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


# Find the AIC-optimal bandwidth and define a search grid around it where we'll fit the 
sim = read.table("~/git/gwr/scratch/sim.txt")
hh = seq(0.13, 0.28, len=S)

indx = index[[i]]
data = sim[indx,]

h = hh[j]
print(h)
model = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), fit.loc=grid, longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005, jacknife=FALSE)
df = sum(sapply(1:length(model[['model']]), function(k) model[['model']][[k]][['tunelist']][['df']]*sum(indx==k) / model[['model']][[k]][['weightsum']]))

#Write LAGR coefficients:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))

#Write LAGR results:
fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, sim[k,2:5])))
dev.resids = gaussian()$dev.resids(sim$Y, fitted, rep(1,N**2))
ll = gaussian()$aic(sim$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))
aic = ll + 2*df
print(aic)

write(c(h, aic), "~/git/gwr/scratch/trace2.txt", append=TRUE)
write.table(coefs, paste("~/git/gwr/scratch/boot/coefs", i, j, "txt", sep="."), row.names=FALSE)


