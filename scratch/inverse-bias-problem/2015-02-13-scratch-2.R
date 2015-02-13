gr = function(x) {
    g = rep(0, 102)
    c = as.matrix(x[1:100])
    d = as.matrix(x[101:102])
    g[1:100] = (t(Q1)%*%T1 + t(Q2)%*%T2)%*%d + (t(Q1)%*%Q1 + t(Q2)%*%Q2 + lambda*Q0)%*%c - t(Q1)%*%obs - t(Q2)%*%obs.
    g[101:102] = (t(T1)%*%T1 + t(T2)%*%T2)%*%d + (t(T1)%*%Q1 + t(T2)%*%Q2)%*%c - t(T1)%*%obs - t(T2)%*%obs.
    
    g = as.vector(2*g)
    g
}

heq = function(x) {
    c = as.matrix(x[1:100])
    t(F1) %*% c
}

heq.jac = function(x) {
    c = as.matrix(x[1:100])
    cbind(t(F1), 0, 0)
}

fn = function(x) {
    c = as.matrix(x[1:100])
    d = as.matrix(x[101:102])
    
    sum((obs - Q1%*%c - T1%*%d)^2) + sum((obs. - Q2%*%c - T2%*%d)^2) + lambda*t(c)%*%Q0%*%c
}

p0 = c(c,d) 
#p0 = rep(0,102)
    
ans = auglag(par=par+rnorm(102,sd=0.5), fn=fn, gr=gr, heq=heq, heq.jac=heq.jac)


c = ans$par[1:100]
d = ans$par[101:102]
f0.smooth.hat = (Q %*% c + T %*% d) 
f0.hat = (Q0 %*% c + T0 %*% d) 

plot(f0.smooth.hat, type='l', bty='n', x=tt, ylim=range(f0(xx)), lwd=2)
par(new=TRUE)
#plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
plot(obs, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0.hat, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0(xx), x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=2, ylim=range(f0(xx)), type='l', lwd=2)

