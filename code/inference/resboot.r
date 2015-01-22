library(lagr)
library(RandomFields)
library(brooks)
library(dplyr)
library(doMC)

registerDoMC(8)

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
set.seed(111882)

#Calculate the coefficient surfaces
B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]

#Seed the RNG
set.seed(2)
indexes = list()
D = list()
epsilon = list()

for (b in 1:5) {

    #The data to use for this simulation setting:
    covariate.model = RMexp(var=1, scale=0.1) + RMnugget(var=0.2)
    d1 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d2 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d3 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]
    d4 = RFsimulate(covariate.model, x=coord, y=coord)@data[[1]]

    D[[b]] = as.matrix(cbind(d1, d2, d3, d4))
    
    loc.x = rep(coord, times=N)
    loc.y = rep(coord, each=N)

    #Simulate the noise term
    epsilon[[b]] = rnorm(N**2, mean=0, sd=1)

    indexes[[b]] = list(sample(400, 100), sample(400, 200), 1:400)
}

process=3 #for (process in 1:3) {
    print(process)

    #Simulation parameters are based on the value of process
    parameters = params[process,]

    #Use the Cholesky decomposition to correlate the random fields:
    S = matrix(parameters[['rho']], 4, 4)
    diag(S) = rep(1, 4)
    L = chol(S)

i=3 #    for (i in 1:3) {
b=3 #        for (b in 1:3) {
            #Force correlation on the Gaussian random fields:
            X = D[[b]] %*% L
            X1 = X[,1]
            X2 = X[,2]
            X3 = X[,3]
            X4 = X[,4]
    
            #Generate the response variable and set up the data.frame:
            eta = X1*B1 + X2*B2 + X3*B3
            Y = eta + epsilon[[b]] * parameters[['sigma']]
            sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), loc.x, loc.y)


            indx = indexes[[b]][[i]]
            data = sim[indx,]
            size = length(indx)
            h = 1.5*size**(-1/6) - 0.36

            file.num = (b-1)*3*6 + (i-1)*6 + process
            write.table(cbind(x=loc.x[indx], y=loc.y[indx], B1=B1[indx], B2=B2[indx], B3=B3[indx], Y=Y[indx]), file=paste("output/truth", file.num, "csv", sep="."))

            #LAGR:
            bw = lagr.tune(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bwselect.method='AIC', bw.type='knn', kernel=epanechnikov, verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005, tol.bw=0.05)
            model = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, varselect.method='AIC', bw=bw, verbose=TRUE, n.lambda=100, lagr.convergence.tol=0.005)

            #Write LAGR coefficients:
            coefs = t(sapply(model[['model']], function(x) x[['coef']]))
            write.table(cbind(x=loc.x[indx], y=loc.y[indx], coefs), file=paste("output/coefs", "lagr", file.num, "csv", sep="."))

            #Write LAGR results:
            fits = as.vector(sapply(model[['model']], function(x) x[['fitted']]))
            resids = Y[indx] - fits
            result = data.frame(x=loc.x[indx], y=loc.y[indx], fitted=fits, residual=resids)
            write.table(result, file=paste("output/results", "lagr", file.num, "csv", sep="."))

            #Clean up:
            rm(model)
            rm(coefs)
            rm(result)
            gc()


            #GWR:
            allvars = replicate(length(indx)**2, c('(Intercept)', 'X1', 'X2', 'X3', 'X4'), simplify=FALSE)
            gwr = lagr(Y~X1+X2+X3+X4, data=data, family='gaussian', coords=c('loc.x','loc.y'), longlat=FALSE, oracle=allvars, bw=h, kernel=epanechnikov, bw.type='knn', verbose=TRUE)

            #Write the GWR coefficients
            cgwr = t(sapply(gwr[['model']], function(x) x[['coef']]))
            write.table(cbind(x=loc.x[indx], y=loc.y[indx], cgwr), file=paste("output/coefs", "gwr", file.num, "csv", sep="."))

            #Write the GWR results:
            fgwr = as.vector(sapply(gwr[['model']], function(x) x[['fitted']]))
            rgwr = Y[indx] - fgwr
            result = data.frame(x=loc.x[indx], y=loc.y[indx], fitted=fgwr, residual=rgwr)
            write.table(result, file=paste("output/results", "gwr", file.num, "csv", sep="."))

            #Clean up:
            rm(gwr)
            rm(cgwr)
            rm(result)
            gc()
        #}
    #}
#}

#As a final act, print the coefficients that we used in the simulation:
zz = range(c(B1, B2, B3))
pdf("~/Desktop/coefs.pdf", 8, 3)
par(oma=c(1,1,1,1))
par('mar'=c(4,1,1,1))
layout(matrix(1:3,1,3))

matplot(B1 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[1]), ylab=NA, cex.lab=2, xrange=zz)
matplot(B2 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[2]), ylab=NA, cex.lab=2, xrange=zz)
matplot(B3 %>% matrix(20,20), border=NA, show.legend=TRUE, axes=FALSE, xlab=expression(beta[3]), ylab=NA, cex.lab=2, xrange=zz)

dev.off()
