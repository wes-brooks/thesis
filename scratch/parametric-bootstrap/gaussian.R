library(lagr)
library(RandomFields)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(3)

#Establish the simulation parameters
tau = 0
sigma.tau = 0
rho = c(rep(0, 2), rep(0.5, 2), rep(0.9, 2)) #rho is the correlation of the covariates
sigma = rep(c(0.5, 1), 3) #sigma is the variance of the noise term
params = data.frame(rho, sigma)

#Generate the covariates:
N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)

#Seed the RNG
set.seed(11181982)

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]

indexes = list(sample(400, 100), sample(400, 200), 1:400)

#The data to use for this simulation setting:
covariate.model = RMexp(var=1, scale=0.1) + RMnugget(var=0.2)
d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)

#Simulate the noise term
epsilon.model = RMexp(var=1, scale=0.3) + RMnugget(var=0.2)
epsilon = RFsimulate(epsilon.model, x=coord, y=coord)@data[[1]]

process=1

#Simulation parameters are based on the value of process
parameters = params[process,]

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 4, 4)
diag(S) = rep(1, 4)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4)) %*% L
X1 = D[,1]
X2 = D[,2]
X3 = D[,3]
X4 = D[,4]

#Generate the response variable and set up the data.frame:
eta = X1*B1 + X2*B2 + X3*B3
Y = eta + epsilon * parameters[['sigma']]
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), loc.x, loc.y)

i = 2
indx = indexes[[i]]
data = sim[indx,]
size = length(indx)
h = 1.5*size**(-1/6) - 0.36


#LAGR:
bw = lagr.tune(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', kernel=epanechnikov, bw.type='knn', tol.bw=0.05, verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005)
model = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005)

#Write LAGR coefficients:
coefs = t(sapply(model[['model']], function(x) x[['coef']]))

#Write LAGR results:
fits = as.vector(sapply(model[['model']], function(x) x[['fitted']]))
resids = Y[indx] - fits
result = data.frame(x=loc.x[indx], y=loc.y[indx], fitted=fits, residual=resids)

