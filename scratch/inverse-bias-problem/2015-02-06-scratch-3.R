f0.smooth = vector()
f1.smooth = vector()
f2.smooth = vector()

f0..smooth = vector()

for (i in 1:100) {
    h = 0.858
    dist = abs(tt-tt[i])
    wt = lagr:::epanechnikov(dist, h)
    X.loc = as.matrix(tt - tt[i])
    
    f0.smooth = c(f0.smooth, lm(f0(tt)~X.loc, weights=wt)$coef[1])
    f1.smooth = c(f1.smooth, lm(f1(tt)~X.loc, weights=wt)$coef[1])
    f2.smooth = c(f2.smooth, lm(f2(tt)~X.loc, weights=wt)$coef[1])
    
    f0..smooth = c(f0..smooth, lm(f0.(tt)~1, weights=wt)$coef[1])
}