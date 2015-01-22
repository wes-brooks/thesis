library(lagr)
library(RandomFields)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(7)

B = 2000

#Generate the covariates:
N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)
n = N**2

#Seed the RNG
set.seed(111882)

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]

#Simulate the covariates:
covariate.model = RMexp(var=1, scale=0.1) + RMnugget(var=0.2)
d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
D = as.matrix(cbind(d1, d2, d3, d4))

#Induce correpation in the covariates:
S = matrix(0.1, 4, 4)
diag(S) = rep(1, 4)
L = chol(S)
X = D %*% L
X1 = X[,1]
X2 = X[,2]
X3 = X[,3]
X4 = X[,4]

#Simulate the noise term
epsilon = rnorm(N**2, mean=0, sd=0.5)

#Generate the response variable and set up the data.frame:
eta = X1*B1 + X2*B2 + X3*B3
Y = eta + epsilon
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4))
sim$loc.x = rep(coord, times=N)
sim$loc.y = rep(coord, each=N)
sim.boot = sim


#estiamte the bandwidth:
bw = lagr.tune(Y~X1+X2+X3+X4, 
               data=sim, 
               family='gaussian', 
               coords=c('loc.x','loc.y'), 
               longlat=FALSE, 
               varselect.method='AICc', 
               bwselect.method='AICc', 
               bw.type='knn', 
               kernel=epanechnikov, 
               verbose=TRUE, 
               n.lambda=100, 
               lagr.convergence.tol=0.005, 
               tol.bw=0.01)

model = lagr(Y~X1+X2+X3+X4, 
             data=sim, 
             family='gaussian', 
             coords=c('loc.x','loc.y'), 
             longlat=FALSE,
             varselect.method='AICc', 
             bw=bw,
             #bw.type='knn',
             #kernel=epanechnikov,
             verbose=TRUE, 
             n.lambda=100, 
             lagr.convergence.tol=0.005)

#######################################
#Bayesian bootstrap Dirichlet weights:
#######################################
wt = list()
for (b in 1:B) {
    raw = sort(c(0, 1, runif(n-1)))
    wt[[b]] = diff(raw) * n
}


#########################################
#residual bootstrap resamples:
#########################################
Y.star = matrix(0, 0, B)

for (j in 1:n) {
    cat(paste("j is ", j, "\n", sep=""))
    Y.star = rbind(Y.star, model$fits[[j]]$model$adamodel$fitted[as.character(j)] + 
                       sample(model$fits[[j]]$model$adamodel$residuals, B, replace=TRUE))
}



AIC.boot.bayes = vector()
df.boot.bayes = vector()

AIC.boot.resid = vector()
df.boot.resid = vector()


#Compute a bandwidth distribution:
indx = which(bw$trace$loss < (min(bw$trace$loss+20)))
max.lik = as.data.frame(bw$trace[indx,1:2])
m.bw = lm(loss ~ log(bw) + I(log(bw)^2), data=max.lik)
mu = -m.bw$coef[2] / m.bw$coef[3] / 2
sd = sqrt(1 / m.bw$coef[3])
bw.boot = bw.b = exp(rnorm(B, mean=mu, sd=sd))

#Run estimation on the bbootstrap draws:
for (b in 1:B) {   
    print(b) 
    print(bw.boot[b])
    
    ####################
    #Bayesian bootstrap:
    ####################
    m.bayes = lagr(Y~X1+X2+X3+X4, 
               data=sim, 
               weights=wt[b]
               family='gaussian', 
               coords=c('loc.x','loc.y'), 
               longlat=FALSE,
               varselect.method='AICc', 
               bw=bw.boot[b],
               bw.type='knn',
               kernel=epanechnikov,
               verbose=TRUE, 
               n.lambda=100, 
               lagr.convergence.tol=0.005)
    
    #Summarize this bootstrap iteration:
    AIC.boot.bayes = c(AIC.boot.bayes, m.bayes$AIC)
    df.boot.bayes = c(df.boot.bayes, m.bayes$df)
    coefs.boot.bayes = t(sapply(m.bayes$fits, function(x) x$coef))
    write.table(coefs.boot.bayes,
                file=paste("/Users/wesley/sim-out/bayesian-bootstrap-coefs-", b, "csv", sep="."))
    
    
    ####################
    #Residual bootstrap:
    ####################
    sim.boot$Y = Y.star[,b]
    m.resid = lagr(Y~X1+X2+X3+X4, 
                   data=sim.boot, 
                   family='gaussian', 
                   coords=c('loc.x','loc.y'), 
                   longlat=FALSE,
                   varselect.method='AICc', 
                   bw=bw.boot[b],
                   bw.type='knn',
                   kernel=epanechnikov,
                   verbose=TRUE, 
                   n.lambda=100, 
                   lagr.convergence.tol=0.005)
    
    #Summarize this bootstrap iteration:
    AIC.boot.resid = c(AIC.boot.resid, m.resid$AIC)
    df.boot.resid = c(df.boot.resid, m.resid$df)
    coefs.boot.resid = t(sapply(m.resid$fits, function(x) x$coef))
    write.table(coefs.boot.resid,
                file=paste("/Users/wesley/sim-out/residual-bootstrap-coefs-", b, "csv", sep="."))
}