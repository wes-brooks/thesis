##Let's try a two-dimensional example:
library(lagr)
library(RandomFields)
library(dplyr)
library(doMC)

registerDoMC(6)

B = 200

#Generate the covariates:
N = 20 #number of width and length divisions in the domain
coord = seq(0, 1, length.out=N)
x.c = rep(coord, each=N)
y.c = rep(coord, times=N)
location = data.frame(x=x.c, y=y.c)
n = N**2

s = 1
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    s = s + as.numeric(args[1])
}



#Seed the RNG
set.seed(s)

#Calculate the coefficient surfaces
test = function(x,y) 3.75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
    3.75*exp(-(9*x+1)^2/49 - (9*y+1)/10) +
    2.5*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
    exp(-(9*x-4)^2 - (9*y-7)^2) - 0.5

t2 = function(x,y) 3*sort(c(0, 5*((x-3/4)^2+(y-1/4)^2) - 1/2,-5*((x-1/4)^2 + (y-3/4)^2) + 1/2 ))[2]

#B1 = RFsimulate(RMexp(var=10, scale=1), x=coord, y=coord)@data[[1]]
#B2 = RFsimulate(RMexp(var=1, scale=1), x=coord, y=coord)@data[[1]]
#B3 = RFsimulate(RMexp(var=0.1, scale=1), x=coord, y=coord)@data[[1]]
B1 = test(location$x, location$y)
B2 = sapply(1:nrow(location), function(k) t2(location$x[k], location$y[k]))
B3 = rep(-0.1, n)

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

#estimate the bandwidth:
bw = lagr.tune(Y~X1+X2+X3+X4, 
               data=sim, 
               family='gaussian', 
               coords=c('x','y'), 
               longlat=FALSE, 
               varselect.method='AICc', 
               bwselect.method='AICc', 
               bw.type='knn', 
               kernel=epanechnikov, 
               verbose=FALSE, 
               n.lambda=100, 
               lagr.convergence.tol=0.005, 
               tol.bw=0.01,
               range=c(0,0.75))

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
    

#Set up a thin-plate spline basis:
theta = function(m, d) {
    #even:
    if (d %% 2 == 0) {
        res = (-1)^(d/2 + m + 1) / 2^(2*m-1) / pi^(d/2) / factorial(m-1) / factorial(m - d/2)
    }
    
    #odd:
    if (d %% 2 ==1) {
        res = gamma(d/2 - m) / 2^(2*m) / pi^(d/2) / factorial(m-1)
    }
    
    res
}


eta = function(r, m, d) {
    #even:
    if ((d %% 2)==0) {
        res = theta(m,d) * r^(2*m-d) * log(r)
    }
    
    #odd:
    if ((d %% 2)==1) {
        res = theta(m,d) * r^(2*m-d)
    }
    
    res
}


W1 = matrix(0, n, n)
for (i in 1:n) {
    h = model$fits[[i]]$bw
    w = epanechnikov(sqrt((location$x[i]-location$x)**2 + (location$y[i]-location$y)**2), h)
    X = as.matrix(cbind(1, location$x-location$x[i], location$y-location$y[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W1[i,] = w * XtX.I[1,] %*% t(X)
}

Q0 = matrix(0, n, n)
for (i in 1:n) {
    Q0[i,] = eta(sqrt((location$x[i]-location$x)**2 + (location$y[i]-location$y)**2), m=2, d=2)
}

T0 = matrix(0, n, 3)
for (i in 1:n) {
    T0[i,1] = 1
    T0[i,2] = location$x[i]
    T0[i,3] = location$y[i]
}

Q1 = W1 %*% Q0
T1 = W1 %*% T0



W2 = matrix(0, n, n)
for (i in 1:n) {
    h = model$fits[[i]]$bw
    w = epanechnikov(sqrt((location$x[i]-location$x)**2 + (location$y[i]-location$y)**2), h)
    X = as.matrix(cbind(1, location$x-location$x[i], location$y-location$y[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W2[i,] = w* XtX.I[2,] %*% t(X)
}

Q2 = W2 %*% Q0
T2 = W2 %*% T0

W3 = matrix(0, n, n)
for (i in 1:n) {
    h = model$fits[[i]]$bw
    w = epanechnikov(sqrt((location$x[i]-location$x)**2 + (location$y[i]-location$y)**2), h)
    X = as.matrix(cbind(1, location$x-location$x[i], location$y-location$y[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W3[i,] = w * XtX.I[3,] %*% t(X)
}

Q3 = W3 %*% Q0
T3 = W3 %*% T0

F1 = (qr(T0) %>% qr.Q(complete=TRUE))[,1:3]
F2 = (qr(T0) %>% qr.Q(complete=TRUE))[,4:400]
R = qr(T0) %>% qr.R

obs = model$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,2) %>% as.matrix
obs.1 = model$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,8) %>% as.matrix
obs.2 = model$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,9) %>% as.matrix

W1.t.I = solve(t(W1))

F2 %*% solve(t(F2) %*% (Q1 + lambda*W1.t.I) %*% F2) %*% t(F2) %*% obs -> c
solve(R) %*% t(F1) %*% (obs - (Q1 - lambda* W1.t.I) %*% c) -> d

f0.smooth.hat = (Q1 %*% c + T1 %*% d) 
f0.hat = (Q0 %*% c + T0 %*% d) 







dir.der = function(x, p) {
    grad = gr(p)[1:400]
    t(x) %*% grad
}

dir.der2 = function(x, p) {
    grad2 = gr2(p[1:400])
    t(x) %*% grad2 %*% x
}    

gr = function(x) {
    g = rep(0, 403)
    c = as.matrix(x[1:400])
    d = as.matrix(x[401:403])
    g[1:400] = (t(Q1)%*%T1 + t(Q2)%*%T2 + t(Q3)%*%T3)%*%d + (t(Q1)%*%Q1 + t(Q2)%*%Q2 + t(Q3)%*%Q3 + lambda*Q0)%*%c - t(Q1)%*%obs - t(Q2)%*%obs.1 - t(Q3)%*%obs.2
    g[401:403] = (t(T1)%*%T1 + t(T2)%*%T2 + t(T3)%*%T3)%*%d + (t(T1)%*%Q1 + t(T2)%*%Q2 + t(T3)%*%Q3)%*%c - t(T1)%*%obs - t(T2)%*%obs.1 - t(T3)%*%obs.2
    
    g = as.vector(2*g)
    g
}

gr2 = function(x) {
    t(Q1)%*%Q1 + t(Q2)%*%Q2 + t(Q3)%*%Q3 + lambda*Q0
}

heq = function(x) {
    c = as.matrix(x[1:400])
    t(F1) %*% c
}

heq.jac = function(x) {
    c = as.matrix(x[1:400])
    cbind(t(F1), 0, 0)
}

fn = function(x) {
    c = as.matrix(x[1:400])
    d = as.matrix(x[401:403])
    
    sum((obs - Q1%*%c - T1%*%d)^2) + sum((obs.1 - Q2%*%c - T2%*%d)^2) + sum((obs.2 - Q3%*%c - T3%*%d)^2) + lambda*t(c)%*%Q0%*%c
}