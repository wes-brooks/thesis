#\beta_0 plots:

plot(x=tt, y=f0(tt), type='l', xlab='t', ylab='response', ylim=bb0, bty='n')
par(new=TRUE)
plot(x=tt, y=f0(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f0.second.derivative(tt), type='l', ylim=bb0, col='red')
par(new=TRUE)
plot(x=tt, y=2*sqrt(3/5) + f0(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f0.second.derivative(tt), type='l', ylim=bb0, col='red', lty=2)
par(new=TRUE)
plot(x=tt, y=-2*sqrt(3/5) + f0(tt)+sapply(m[['model']], function(x) x[['bw']]) / 10 * f0.second.derivative(tt), type='l', ylim=bb0, col='red', lty=2)


par(new=TRUE)
plot(x=tt, y=coefs[,1], type='l', xlab='t', ylab='response', ylim=bb0, bty='n', col='blue')

sdev = sqrt(diag(SS)[(0:99)*3+1])
par(new=TRUE)
plot(x=tt, y=2*sdev + coefs[,1], type='l', ylim=bb0, col='blue', lty=2)
par(new=TRUE)
plot(x=tt, y=-2*sdev + coefs[,1], type='l', ylim=bb0, col='blue', lty=2)

par(new=TRUE)
plot(x=tt, y=beta.star.4[[2]][2,], type='l', xlab='t', ylab='response', ylim=bb1, bty='n', col='green')


par(new=TRUE)
plot(x=tt, y=beta.star.2[[10]][2,], type='l', xlab='t', ylab='response', ylim=bb1, bty='n', col='yellow', lwd=2)

