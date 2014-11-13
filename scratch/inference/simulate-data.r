library(RandomFields)

#Establish the simulation parameters
tau = 0
sigma.tau = 0
rho = 0.5
sigma = 0.5

#Generate the covariates:
N = 14
coord = seq(0, 1, length.out=N)
loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)
grid = cbind(loc.x, loc.y)

set.seed(123)

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]

#The data to use for this simulation setting:
covariate.model = RMexp(var=1, scale=0.1) + RMnugget(var=0.2)
d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
D = cbind(d1, d2, d3, d4)

#Simulate random error
epsilon = rnorm(N**2, mean=0, sd=1)

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(rho, 4, 4)
diag(S) = rep(1, 4)
L = chol(S)

#Force correlation on the Gaussian random fields:
X = D %*% L
X1 = X[,1]
X2 = X[,2]
X3 = X[,3]
X4 = X[,4]

#Generate the response variable and set up the data.frame:
eta = X1*B1 + X2*B2 + X3*B3
Y = eta + epsilon * sigma
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), loc.x, loc.y)
