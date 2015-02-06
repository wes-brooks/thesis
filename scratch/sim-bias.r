library(lagr)
library(RandomFields)
library(dplyr)
library(doMC)

registerDoMC(6)

#Generate the covariates:
N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)
x.c = rep(coord, each=N)
y.c = rep(coord, times=N)
location = data.frame(x=x.c, y=y.c)
n = N**2
coefs = array(dim=c(400,5,200))
bw=0.08


#Calculate the coefficient surfaces
test = function(x,y) 3.75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
    3.75*exp(-(9*x+1)^2/49 - (9*y+1)/10) +
    2.5*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
    exp(-(9*x-4)^2 - (9*y-7)^2) - 0.5

t2 = function(x,y) 3*sort(c(0, 5*((x-3/4)^2+(y-1/4)^2) - 1/2,-5*((x-1/4)^2 + (y-3/4)^2) + 1/2 ))[2]

B2.. = function(x,y) {
    f.id = apply(cbind(x,y), 1, function(x) order(c(0, 5*((x[1]-3/4)^2+(x[2]-1/4)^2) - 1/2,-5*((x[1]-1/4)^2 + (x[2]-3/4)^2) + 1/2 ))[2])
    ifelse(f.id == 1, 0, ifelse(f.id==2, 30, -30))
}

B1 = test(location$x, location$y)
B2 = sapply(1:nrow(location), function(k) t2(location$x[k], location$y[k]))
B3 = rep(-0.1, n)


for (s in 1:200) {
    print(paste("beginning", s))
    set.seed(s)

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
    epsilon = rnorm(N**2, mean=0, sd=1)
    
    #Generate the response variable and set up the data.frame:
    eta = X1*B1 + X2*B2 + X3*B3
    Y = eta + epsilon
    sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4))
    sim$x = x.c
    sim$y = y.c
    
    model = lagr(Y~X1+X2+X3+X4, 
                 data=sim, 
                 family='gaussian', 
                 coords=c('x','y'), 
                 longlat=FALSE,
                 varselect.method='AICc', 
                 bw=bw,
                 bw.type='knn',
                 kernel=epanechnikov,
                 verbose=FALSE, 
                 n.lambda=100, 
                 lagr.convergence.tol=0.005)
 
    
    cc = model$fits %>% sapply(function(x) if(!is.null(x$coef)) x$coef else rep(NA,5)) %>% t
    coefs[,,s] = cc
    
    #Empirical second derivative of beta.0
    slx = sapply(model$fits, function(x) x$model$adamodel$coef[8])
    sly = sapply(model$fits, function(x) x$model$adamodel$coef[9])
    cc.0 = vector()
    for (i in 1:N**2) {
        h = model$fits[[i]]$bw
        X = cbind(location$x - location$x[i], location$y - location$y[i])
        w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
        
        m2x = lm(slx~X, weights=w)
        m2y = lm(sly~X, weights=w)
        cc.0 = c(cc.0, m2x$coef[2] + m2y$coef[3])
    }
    empirical.bias.0 = cc.0 * sapply(model$fits, function(x) x$bw**2) / 10
    
    cc.1 = vector()
    for (i in 1:N**2) {
        h = model$fits[[i]]$bw / 2
        X = cbind(location$x - location$x[i], location$y - location$y[i])
        w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
        
        m2x = lm(slx~X, weights=w)
        m2y = lm(sly~X, weights=w)
        cc.1 = c(cc.1, m2x$coef[2] + m2y$coef[3])
    }
    empirical.bias.1 = cc.1 * sapply(model$fits, function(x) x$bw**2) / 10
    
    
    cc.2 = vector()
    for (i in 1:N**2) {
        h = model$fits[[i]]$bw / 4
        X = cbind(location$x - location$x[i], location$y - location$y[i])
        w = lagr:::epanechnikov(sqrt(apply(X, 1, function(x) sum(x**2))), h)
        
        m2x = lm(slx~X, weights=w)
        m2y = lm(sly~X, weights=w)
        cc.2 = c(cc.2, m2x$coef[2] + m2y$coef[3])
    }
    empirical.bias.2 = cc.2 * sapply(model$fits, function(x) x$bw**2) / 10
    
    pdf(paste("~/Desktop/bias-plots/bias.", s, ".pdf", sep=""), width=18,height=6)
    empirical.bias.0 %>% plot(type='l', col='red', ylim=c(-2,2))
    par(new=TRUE)
    empirical.bias.1 %>% plot(type='l', col='blue', ylim=c(-2,2))
    par(new=TRUE)
    empirical.bias.2 %>% plot(type='l', col='green', ylim=c(-2,2))
    par(new=TRUE)
    (cc[,2]-B1) %>% plot(type='l', ylim=c(-2,2))
    dev.off()
}




#Smooth the expected bias
n=400
squaredist = matrix(0,n,n)
D.coords = as.matrix(cbind(location$x, location$y))
q = 2
for (c in 1:q) {
    this.coord = matrix(rep(D.coords[,c], times=n), n, n)
    squaredist = squaredist + (this.coord-t(this.coord))**2
}
D = sqrt(squaredist)

bias.1 = model$fits %>% sapply(function(x) x$bw^2) * B1..(location$x, location$y) / 10
b.smoothed.1 = vector()
B1.smooth = vector()
for (i in 1:400) {
    h = model$fits[[i]]$bw
    dist = D[i,]
    wt = lagr:::epanechnikov(dist, h)
    X.loc = cbind(location$x-location$x[i], location$y-location$y[i])
    
    b.smoothed.1 = c(b.smoothed.1, sum(wt*bias.1) / sum(wt))
    B1.smooth = c(B1.smooth, lm(B1~X.loc, weights=wt)$coef[1])
}


#Plot the observed bias versus the theoretical bias:
(model$fits %>% sapply(function(x) x$bw^2) * B1..(location$x, location$y) / 10 ) %>% plot(type='l', ylim=c(-2,2))
par(new=TRUE)
(sapply(1:400, function(i) mean(coefs[i,2,], na.rm=TRUE))-B1) %>% plot(type='l', col='red', ylim=c(-2,2))
par(new=TRUE)
b.smoothed.1 %>% plot(type='l', col='blue', ylim=c(-2,2))
legend(c("Theoretical asymptotic bias", "mean bias in 80 realizations", "Smoothed theoretical bias"), col=c("black", "red", 'blue'), x="bottomright", bty='n', lty=c(1,1,1))















#Smooth the expected bias
n=400
squaredist = matrix(0,n,n)
D.coords = as.matrix(cbind(location$x, location$y))
q = 2
for (c in 1:q) {
    this.coord = matrix(rep(D.coords[,c], times=n), n, n)
    squaredist = squaredist + (this.coord-t(this.coord))**2
}
D = sqrt(squaredist)

bias.2 = model$fits %>% sapply(function(x) x$bw^2) * B2..(location$x, location$y) / 10
b.smoothed.2 = vector()
B2.smooth = vector()
for (i in 1:400) {
    h = model$fits[[i]]$bw
    dist = D[i,]
    wt = lagr:::epanechnikov(dist, h)
    X.loc = cbind(location$x[i]-location$x, location$y[i]-location$y)
    
    b.smoothed.2 = c(b.smoothed.2, sum(wt*bias.2) / sum(wt))
    B2.smooth = c(B2.smooth, lm(B2~X.loc, weights=wt)$coef[1])
}


#Plot the observed bias versus the theoretical bias:
(model$fits %>% sapply(function(x) x$bw^2) * B2..(location$x, location$y) / 10 ) %>% plot(type='l', ylim=c(-2,2))
par(new=TRUE)
(sapply(1:400, function(i) mean(coefs[i,3,], na.rm=TRUE))-B2.smooth) %>% plot(type='l', col='red', ylim=c(-2,2))
par(new=TRUE)
b.smoothed.2 %>% plot(type='l', col='blue', ylim=c(-2,2))
legend(c("Theoretical asymptotic bias", "mean bias in 80 realizations", "Smoothed theoretical bias"), col=c("black", "red", 'blue'), x="bottomright", bty='n', lty=c(1,1,1))

