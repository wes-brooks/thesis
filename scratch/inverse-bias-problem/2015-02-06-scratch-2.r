hh = seq(0.8,0.9, len=6)
co = list()
co.sh = list()
for (j in 1:length(hh)) {
    h = hh[j]
    print(h)
    m = lagr(y~x1+x2, data=df, family='gaussian', coords='t', varselect.method='wAICc', bw=h, kernel=epanechnikov, bw.type='dist', verbose=TRUE, lagr.convergence.tol=0.005, lambda.min.ratio=0.01, n.lambda=80)

    co.sh[[j]] = sapply(m$fits, function(x) x$coef) %>% t
    co[[j]] = sapply(m$fits, function(x) x$model$adamodel$coef) %>%t
}

proj1 = vector()
proj2 = vector()
for (i in 1:100) {
    df1 = data.frame(x=hh, y=sapply(co, function(x) x[i,2]))
    df2 = data.frame(x=hh, y=sapply(co.sh, function(x) x[i,2]))
        
    proj1 = c(proj1, lm(y~I(x^2), data=df1)$coef[1])
    proj2 = c(proj2, lm(y~I(x^2), data=df2)$coef[1])
}

range(proj1, proj2, f1(tt)) -> yy
plot(x=tt, y=f1(tt), type='l', col='red', ylim=yy, bty='n', xlab="t", ylab="B2")
par(new=TRUE)
plot(x=tt, y=proj1, type='l', col='black', lty=2, ylim=yy, bty='n', xaxt='n', yaxt='n', ann=FALSE)
par(new=TRUE)
plot(x=tt, y=proj2, type='l', col='black', lty=3, ylim=yy, bty='n', xaxt='n', yaxt='n', ann=FALSE)
par(new=TRUE)
plot(x=tt, y=co.sh[[3]][,2], type='l', col='blue', ylim=yy, bty='n', xaxt='n', yaxt='n', ann=FALSE)
