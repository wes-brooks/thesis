b = bw[['trace']][order(bw[['trace']][,1]),c(1,2)]
spline = smooth.spline(b)

xxx = seq(b[1,1], tail(b[,1],1), length.out=1001)
smooth = predict(spline, xxx)

maxi = which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))>0.9999)[1]
mini = tail(which(cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))<0.0001),1)

xxx = seq(xxx[mini], xxx[maxi], length.out=1001)
smooth = predict(spline, xxx)

pp = cumsum(exp(-0.5*smooth$y))/sum(exp(-0.5*smooth$y))

bws = xxx[sapply(runif(19), function(x) which(x<pp)[1])]
models = list()
models[[1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], N=100, longlat=FALSE, varselect.method='AICc', bw=bw.lagr[['bw']], kernel=epanechnikov, bw.type='knn', simulation=TRUE, verbose=TRUE)

for (bw in bws) {
    models[[length(models)+1]] = lagr(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', coords=sim[,c('loc.x','loc.y')], N=100, longlat=FALSE, varselect.method='AICc', bw=bw, kernel=epanechnikov, bw.type='knn', simulation=TRUE, verbose=TRUE)
}
