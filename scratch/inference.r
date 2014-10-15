library(lagr)
# library(RandomFields)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(8)

#Establish the simulation parameters
# tau = 0
# sigma.tau = 0
# rho = c(rep(0, 2), rep(0.5, 2), rep(0.9, 2)) #rho is the correlation of the covariates
# sigma = rep(c(0.5, 1), 3) #sigma is the variance of the noise term
# params = data.frame(rho, sigma)

#Generate the covariates:
N = 14 # number of width and length divisions in the domain
B = 200 # number of bootstrap resamplings of the original data
S = 80 # steps in the grid search for bandwidth
coord = seq(0, 1, length.out=N)
# 
# #Seed the RNG
set.seed(11181982)
# 
# #Calculate the coefficient surfaces
# B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
# B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
# B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]
# 
# #The data to use for this simulation setting:
# covariate.model = RMexp(var=1, scale=0.1) + RMnugget(var=0.2)
# d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
# d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
# d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
# d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]


# #Simulate the noise term
# epsilon.model = RMexp(var=1, scale=0.5) + RMnugget(var=0.2)
# epsilon = RFsimulate(epsilon.model, x=coord, y=coord)@data[[1]] #= rnorm(N**2, mean=0, sd=1)

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



# 
# process = 1    
#     #Simulation parameters are based on the value of process
#     parameters = params[process,]
#     
#     #Use the Cholesky decomposition to correlate the random fields:
#     Sigma = matrix(parameters[['rho']], 4, 4)
#     diag(Sigma) = rep(1, 4)
#     L = chol(Sigma)
#     
#     #Force correlation on the Gaussian random fields:
#     D = as.matrix(cbind(d1, d2, d3, d4)) %*% L
#     X1 = D[,1]
#     X2 = D[,2]
#     X3 = D[,3]
#     X4 = D[,4]
#     
#     #Generate the response variable and set up the data.frame:
#     eta = X1*B1 + X2*B2 + X3*B3
#     Y = eta + epsilon * parameters[['sigma']]
#     sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), loc.x, loc.y)

    # Find the AIC-optimal bandwidth and define a search grid around it where we'll fit the 
    #bw = lagr.tune(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', kernel=epanechnikov, bw.type='knn', tol.bw=0.01, verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005)
sim = read.table("~/git/gwr/scratch/sim.txt")
hh = seq(0.1, 0.5, len=S)
    
        indx = index[[i]]
        data = sim[indx,]
        
            h = hh[j]
            print(h)
            model = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005, jacknife=TRUE)
            df = sum(sapply(model[['model']], function(x) tail(x[['tunelist']][['df-local']], 1)))
            
            #Write LAGR coefficients:
            coefs = t(sapply(model[['model']], function(x) x[['coef']]))
            
            #Write LAGR results:
            fitted = sapply(1:196, function(k) sum(coefs[k,] * cbind(1, data[k,2:5])))
            dev.resids = gaussian()$dev.resids(data$Y, fitted, rep(1,N**2))
            ll = gaussian()$aic(data$Y, N**2, fitted, rep(1,N**2), sum(dev.resids))
            err = ll
            #aic[[i]] = c(aic[[i]], ll + 2*df)
            print(ll)

    write(c(h, err), "~/git/gwr/scratch/trace.txt", append=TRUE)

