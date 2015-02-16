tt = seq(0,5,len=100)

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


eta. = function(r, m, d) {
    #even:
    if ((d %% 2) == 0) {
        res = theta(m,d) * ((2*m-d)*r^(2*m-d-1)*log(r) + r^(2*m-d-1))
    }
    
    #odd:
    if ((d %% 2) == 1) {
        res = theta(m,d) * (2*m-d)*r^(2*m-d-1)
    }
    
    res
}



W1 = matrix(0, n, n)
for (i in 1:n) {
    #h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = as.matrix(cbind(1, tt-tt[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W1[i,] = w * XtX.I[1,] %*% t(X)
}

Q0 = matrix(0, n, n)
for (i in 1:n) {
    Q0[i,] = eta(abs(tt-tt[i]), m=2, d=1)
}

T0 = matrix(0, n, 2)
for (i in 1:n) {
    T0[i,1] = 1
    T0[i,2] = tt[i]
}

Q1 = W1 %*% Q0
T1 = W1 %*% T0



W2 = matrix(0, n, n)
for (i in 1:n) {
    #h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = as.matrix(cbind(1, tt-tt[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W2[i,] = w* XtX.I[2,] %*% t(X)
}

# Q0. = matrix(0, n, n)
# for (i in 1:n) {
#     Q0.[i,] = eta.(abs(tt-tt[i]), m=2, d=1)
# }
# 
# T0. = matrix(0, n, 2)
# for (i in 1:n) {
#     T0.[i,1] = 0
#     T0.[i,2] = 1
# }

Q2 = W2 %*% Q0
T2 = W2 %*% T0


xx = seq(0,5,len=100)
F1 = (qr(T0) %>% qr.Q(complete=TRUE))[,1:2]
F2 = (qr(T0) %>% qr.Q(complete=TRUE))[,3:100]
R = qr(T0) %>% qr.R

obs = m$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,1) %>% as.matrix
obs. = m$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,4) %>% as.matrix

QtI = solve(t(W1)) %*% solve(t(Q0))
#Iterative solution for Y, Z:
d = matrix(0,2,1)
F2 %*% solve(t(F2)%*%(t(W)%*%W + t(W.)%*%W.) %*%Q0%*% F2 + diag(rep(lambda, 98))) %*% t(F2) %*% (W%*%obs + W.%*%obs. - (t(W)%*%W + t(W.)%*%W.)%*%T0%*%d) -> c
solve(t(F1) %*% (t(W)%*%W + t(W.)%*%W.) %*% F1 %*% R) %*% t(F1) %*% (t(W)%*%obs + t(W.)%*%obs. - (t(W)%*%W + t(W.)%*%W.)%*%Q0%*%c) -> d
d

lambda=0.01
#F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% ginv(Q) %*% Q0 %*% F2) %*% t(F2) %*% f0.smooth -> c
#solve(R) %*% (t(F1) %*% f0.smooth - t(F1) %*% Q %*% c - lambda*t(F1) %*% ginv(Q) %*% Q0 %*%c) -> d

#F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% solve(W) %*% F2) %*% t(F2) %*% f0.smooth -> c
#solve(R) %*% (t(F1) %*% f0.smooth - t(F1) %*% Q %*% c - lambda*t(F1) %*% solve(W) %*%c) -> d
gcv = vector()
for (lambda in exp(seq(log(0.1), log(0.00001), len=200))) {
    #F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% solve(W) %*% F2) %*% t(F2) %*% obs -> c
    #solve(R) %*% (t(F1) %*% obs - t(F1) %*% Q %*% c - lambda*t(F1) %*% solve(W) %*% c) -> d
    F2 %*% solve(t(F2) %*% (Q + lambda*QtI%*%Q0) %*% F2) %*% t(F2) %*% obs -> c
    solve(R) %*% t(F1) %*% (obs - (Q - lambda* QtI%*%Q0) %*% c) -> d
    
    
    f0.smooth.hat = (Q1 %*% c + T1 %*% d) 
    f0.hat = (Q0 %*% c + T0 %*% d) 
    f0.slope.hat.smooth = (Q2 %*% c + T2 %*% d) 
    
    #GCV criterion:
    A.c = F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% solve(W) %*% F2) %*% t(F2)
    A.d = solve(R) %*% (t(F1) - t(F1) %*% Q %*% A.c - lambda*t(F1) %*% solve(W) %*% A.c)
    A = Q %*% A.c + T %*% A.d
    gcv = c(gcv, sum((obs-f0.smooth.hat)^2) / (n-sum(diag(A)))^2)
}


f0.smooth.hat = (Q1 %*% c + T1 %*% d) 
f0.hat = (Q0 %*% c + T0 %*% d) 
f0.slope.hat.smooth = (Q2 %*% c + T2 %*% d) 


plot(f0.smooth.hat, type='l', bty='n', x=tt, ylim=range(f0(xx)), lwd=2)
par(new=TRUE)
#plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
plot(obs, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0.hat, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0(xx), x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=2, ylim=range(f0(xx)), type='l', lwd=2)



plot(f0.slope.hat.smooth, type='l', bty='n', x=tt, ylim=c(-15,5), lwd=2)
par(new=TRUE)
#plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
plot(obs., x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=c(-15,5), type='l', lwd=2)
par(new=TRUE)
plot(f0..smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=c(-15,5), type='l', lwd=2)

