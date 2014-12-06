//\beta_1 plots

plot(x=tt, y=f1(tt), type='l', xlab='t', ylab='response', ylim=bb1, bty='n')
par(new=TRUE)
plot(x=tt, y=f1(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f1.second.derivative(tt), type='l', ylim=bb1, col='red')
par(new=TRUE)
plot(x=tt, y=2*sqrt(3/5/11) + f1(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f1.second.derivative(tt), type='l', ylim=bb1, col='red', lty=2)
par(new=TRUE)
plot(x=tt, y=-2*sqrt(3/5/11) + f1(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f1.second.derivative(tt), type='l', ylim=bb1, col='red', lty=2)


par(new=TRUE)
plot(x=tt, y=coefs[,2], type='l', xlab='t', ylab='response', ylim=bb1, bty='n', col='blue')

sdev = sqrt(diag(SS)[(0:99)*3+2])
par(new=TRUE)
plot(x=tt, y=2*sdev + coefs[,2], type='l', ylim=bb1, col='blue', lty=2)
par(new=TRUE)
plot(x=tt, y=-2*sdev + coefs[,2], type='l', ylim=bb1, col='blue', lty=2)

par(new=TRUE)
plot(apply(sapply(beta.star.4, function(x) x[2,]), 1, function(y) quantile(y, 0.2)), type='l', col='green', ylim=bb1)
par(new=TRUE)
plot(apply(sapply(beta.star.4, function(x) x[2,]), 1, function(y) quantile(y, 0.8)), type='l', col='green', ylim=bb1)


par(new=TRUE)
plot(x=tt, y=beta.star.4[[2]][2,], type='l', xlab='t', ylab='response', ylim=bb1, bty='n', col='green')


par(new=TRUE)
plot(x=tt, y=beta.star.2[[10]][2,], type='l', xlab='t', ylab='response', ylim=bb1, bty='n', col='yellow', lwd=2)


