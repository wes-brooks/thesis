library(sp, lib.loc"R-libs")
library(iterators, lib.loc="R-libs")
library(foreach, lib.loc="R-libs")
library(pryr, lib.loc="R-libs")
library(Matrix)
library(lagr, lib.loc="R-libs")

#Generate the covariates:
N = 14 # number of width and length divisions in the domain
B = 2 # number of bootstrap resamplings of the original data
S = 5 # steps in the grid search for bandwidth
hh = seq(0.13, 0.4, len=S)
coord = seq(0, 1, length.out=N)

# #Seed the RNG
set.seed(11181982)

#Generate the bootstrap draws:
index = list()
index[[1]] = 1:(N^2-1)
if (B > 1) {
    for (i in 2:B) {
        index[[i]] = sample(1:(N^2-1), replace=TRUE)
    }
}

if (interactive()) {
    i=1
    j=1
} else {
    iter = commandArgs(trailingOnly=TRUE)[1] 
    iter = as.numeric(iter)
    i = ((iter-1) %/% S) + 1
    j = ((iter-1) %% S) + 1

    cat(paste("Iteration: ", iter, "; bootstrap: ", i, ", bandwidth: ", hh[j], "\n", sep=""))
    
    #On the first pass, write the indicies
    if (iter==1) {
        ind = as.data.frame(index)
        colnames(ind) = NULL
        write.table(ind, "output/indexes.txt", row.names=FALSE, col.names=FALSE)
    }
}

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
grid = cbind(loc.x, loc.y)

sim = read.table("sim.txt")

h = hh[j]
model = lagr(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005, jacknife=TRUE, bootstrap.index=index[[i]])

#Write jacknife fitting results:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))
fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, sim[k,2:5])))
dev.resids = gaussian()$dev.resids(sim$Y, fitted, rep(1,N**2))
ll = gaussian()$aic(sim$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))

write(c(h, ll), paste("output/trace", "jacknife", i, j, "txt", sep="."), append=TRUE)
cat(paste("Bandwith: ", h, "; Jacknife loss: ", ll, "\n", sep=''))

#Write jacknife LAGR coefficients:
for (k in 1:(N**2)) {
    beta = model[['model']][[k]][['model']][['beta']][1:5,]
    models = as.factor(apply(beta, 2, function(x) paste(as.integer(x!=0), collapse="")))
    aic.trace = model[['model']][[k]][['model']][['results']][['AIC']]
    
    best = vector()
    for (l in levels(models)) {
        best = c(best, which.min(ifelse(models==l, aic.trace, Inf)))
    }
    
    write.table(beta[,best], paste("output/beta", "jacknife", i, j, k, "txt", sep="."), row.names=FALSE, col.names=FALSE)
    write(aic.trace[best], paste("output/aic", "jacknife", i, j, k, "txt", sep="."))
}


#Anti-jacknife:
model = lagr(Y~X1+X2+X3+X4, data=sim, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=FALSE, n.lambda=100, lagr.convergence.tol=0.005, jacknife='anti', bootstrap.index=index[[i]])

#Write anti-jacknife fitting results:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))
fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, sim[k,2:5])))
dev.resids = gaussian()$dev.resids(sim$Y, fitted, rep(1,N**2))
ll = gaussian()$aic(sim$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))
df = sum(sapply(1:196, function(k) model[['model']][[k]][['tunelist']][['df-local']]))
aic = ll + 2*df

#Use the AIC for anti-jacknife bandwidth tuning:
write(c(h, aic), paste("output/trace", "anti", i, j, "txt", sep="."), append=TRUE)
cat(paste("Bandwith: ", h, "; Anti-jacknife AIC: ", aic, "\n", sep=''))


#Write anti-jacknife LAGR coefficients:
for (k in 1:(N**2)) {
    beta = model[['model']][[k]][['model']][['beta']][1:5,]
    models = as.factor(apply(beta, 2, function(x) paste(as.integer(x!=0), collapse="")))
    aic.trace = model[['model']][[k]][['model']][['results']][['AIC']]
    
    best = vector()
    for (l in levels(models)) {
        best = c(best, which.min(ifelse(models==l, aic.trace, Inf)))
    }
    
    write.table(beta[,best], paste("output/beta", "anti", i, j, k, "txt", sep="."), row.names=FALSE, col.names=FALSE)
    write(aic.trace[best], paste("output/aic", "anti", i, j, k, "txt", sep="."))
}