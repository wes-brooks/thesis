library(lagr)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(8)

#Generate the covariates:
N = 14 # number of width and length divisions in the domain
B = 25 # number of bootstrap resamplings of the original data
S = 61 # steps in the grid search for bandwidth
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
    i=6
    j=10
} else {
    iter = commandArgs(trailingOnly=TRUE)[1] 
    print(iter)
    iter = as.numeric(iter)
    i = ((iter-1) %/% S) + 1
    j = ((iter-1) %% S) + 1
    
    print(i)
    print(j)
    
    #On the first pass, write the indicies
    if (iter==1) {
        ind = as.data.frame(index)
        colnames(ind) = NULL
        write.table(ind, "~/git/gwr/output/indexes.txt", row.names=FALSE, col.names=FALSE)
    }
}

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
grid = cbind(loc.x, loc.y)

# Import the simulated data
sim = read.table("~/git/gwr/scratch/sim.txt")
hh = seq(0.13, 0.43, len=S)

indx = index[[i]]
data = sim[indx,]

h = hh[j]
print(h)
model = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), fit.loc=grid, longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005, jacknife=FALSE)

#Write LAGR coefficients:
for (k in 1:(N**2)) {
    beta = model[['model']][[k]][['model']][['beta']][1:5,]
    models = apply(beta, 2, function(x) paste(as.integer(x!=0), collapse="")) %>% as.factor
    aic.trace = model[['model']][[k]][['model']][['results']][['AIC']]
    
    best = vector()
    for (l in levels(models)) {
        best = c(best, which.min(ifelse(models==l, aic.trace, Inf)))
    }
    
    write.table(beta[,best], paste("~/git/gwr/output/beta", i, j, k, "txt", sep="."))
    write(aic.trace[best], paste("~/git/gwr/output/aic", i, j, k, "txt", sep="."))
}

#Write LAGR results:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))
fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, sim[k,2:5])))
df = sum(sapply(1:length(model[['model']]), function(k) model[['model']][[k]][['tunelist']][['df']]*sum(indx==k)/model[['model']][[k]][['weightsum']]))
dev.resids = gaussian()$dev.resids(data$Y, fitted[indx], rep(1,N**2))
ll = gaussian()$aic(data$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))
aic = ll + 2*df
print(aic)

write(c(h, aic), "~/git/gwr/output/trace.txt", append=TRUE)


