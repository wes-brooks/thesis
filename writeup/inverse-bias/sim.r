library(lagr)
library(doMC)
registerDoMC(7)

## @knitr simulate-data
#Set constants:
B1 = 1
n = 100
tt = seq(0, 5, len=n)

#Coefficient functions:
f0 = function(x) {-0.2*x^4 + x^3 + 0.7*(x-1)^2 - 4*(x-2) + 1}
f1 = function(x) {0.2*x^2 - x + cos(x)}
f2 = function(x) {(x+1)^(-1) - 1/6}

#Generate the raw data sets (not bootstrap)
set.seed(33)
#simulate data
X1 = rnorm(n, mean=3, sd=2)
X2 = rnorm(n, mean=3, sd=2)
y = f0(tt) + X1*f1(tt) + X2*f2(tt) + rnorm(n)
df = data.frame(y, x1=X1, x2=X2, t=tt)
  
#Fit a VCR model to the simulated data by LAGR
m = lagr(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc', bw=0.858, kernel=epanechnikov, bw.type='dist', verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)

#Store the fitted values, residuals, and coefficient estimates:
fitted = sapply(m[['fits']], function(x) tail(x[['fitted']],1))
resid = df$y - fitted
coefs = t(sapply(m$fits, function(x) x[['model']][['beta']][,ncol(x[['model']][['beta']])]))


f0.smooth = vector()
f1.smooth = vector()
f2.smooth = vector()

for (i in 1:100) {
    h = 0.858
    dist = abs(tt-tt[i])
    wt = lagr:::epanechnikov(dist, h)
    X.loc = as.matrix(tt - tt[i])
    
    f0.smooth = c(f0.smooth, lm(f0(tt)~X.loc, weights=wt)$coef[1])
    f1.smooth = c(f1.smooth, lm(f1(tt)~X.loc, weights=wt)$coef[1])
    f2.smooth = c(f2.smooth, lm(f2(tt)~X.loc, weights=wt)$coef[1])

    
}
