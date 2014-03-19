require(devtools)
install_github('SGL', 'wrbrooks')
install_github('gwselect', 'wrbrooks', ref='agLasso')
require(SGL)
require(gwselect)
require(geoR)
require(doMC)
registerDoMC(cores=3)

B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(c(0, 0.1), each=9)
rho = rep(rep(c(0, 0.5, 0.8), each=3), times=2)
sigma.tau = rep(c(0, 0.03, 0.1), times=6)
params = data.frame(tau, rho, sigma.tau)

#This is the gradient function:
B1 = matrix(rep(coord, N), N, N)

#Read command-line parameters
cluster=NA
process=2

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]

#Get two (independent) Gaussian random fields:
d1 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d2 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d4 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d5 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))

loc.x = d1$coords[,1]
loc.y = d1$coords[,2]

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1$data, d2$data, d3$data, d4$data, d5$data)) %*% L
    
#
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)

#if (parameters[['function.type']] == 'step') {B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)}
#if (parameters[['function.type']] == 'gradient') {B1 = matrix(rep(1-coord, N), N, N)}


if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=1)}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['sigma.tau']]))$data}

#
mu = X1*B1
Y = mu + epsilon

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fitloc = cbind(rep(seq(0,1, length.out=N), each=N), rep(seq(0,1, length.out=N), times=N))

#Reduce the sample size
indx = sample(1:900, 100)
fitloc = fitloc[indx,]

vars = as.vector(B1!=0)
oracle = list()
for (i in 1:N**2) { 
    if (vars[i]) {
        oracle[[i]] = c("X1")
    } else { 
        oracle[[i]] = character(0)
    }
}

#Find the optimal bandwidth and use it to generate a model:
bw = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, varselect.method="AICc", range=c(0,1), gweight=bisquare, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", parallel=FALSE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson')
bw=0.827034
ml = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, varselect.method='AICc', gweight=bisquare, bw=bw[['bw']], bw.type='knn', parallel=FALSE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson', simulation=TRUE)
