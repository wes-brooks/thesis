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

phi = list()
phi[[1]] = function(t) 1
phi[[2]] = function(t) t

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



W = matrix(0, n, n)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = as.matrix(cbind(1, tt-tt[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W[i,] = w*XtX.I[1,1] + w*XtX.I[2,1]*(tt-tt[i])
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

Q = W %*% Q0
T = W %*% T0



W. = matrix(0, n, n)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    X = as.matrix(cbind(1, tt-tt[i]))
    XtX.I = solve(t(X) %*% diag(w) %*% X)
    
    W.[i,] = w*XtX.I[2,1]*(tt-tt[i]) + w*XtX.I[2,2]*(tt-tt[i])^2
}

Q0. = matrix(0, n, n)
for (i in 1:n) {
    Q0.[i,] = eta.(abs(tt-tt[i]), m=2, d=1)
}

T0. = matrix(0, n, 2)
for (i in 1:n) {
    T0.[i,1] = 0
    T0.[i,2] = 1
}

Q. = W. %*% Q0.
T. = W. %*% T0.


xx = seq(0,5,len=100)
F1 = (qr(T) %>% qr.Q(complete=TRUE))[,1:2]
F2 = (qr(T) %>% qr.Q(complete=TRUE))[,3:100]
R = qr(T) %>% qr.R

obs = m$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,1) %>% as.matrix
obs. = m$fits %>% sapply(function(x) x$model$adamodel$coef) %>% t %>% `[`(,4) %>% as.matrix

QtI = solve(Q0) %*% solve(W)
#Iterative solution for Y, Z:
d = matrix(0,2,1)
F2 %*% solve(t(F2) %*% (Q + lambda*QtI%*%Q0 + QtI%*%t(Q.)%*%Q.) %*% F2) %*% t(F2) %*% (obs + QtI%*%t(Q.)%*%(obs. - T.%*%d) ) -> c
solve(R + t(F1)%*%QtI%*%t(Q.)%*%T.) %*% t(F1) %*% (obs + Q%*%c + QtI%*%t(Q.)%*%obs. + QtI%*%(t(Q.)%*%Q. + lambda*Q0) %*% c) -> d


lambda=0.1
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
    
    
    f0.smooth.hat = (Q %*% c + T %*% d) 
    f0.hat = (Q0 %*% c + T0 %*% d) 
    
    #GCV criterion:
    A.c = F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% solve(W) %*% F2) %*% t(F2)
    A.d = solve(R) %*% (t(F1) - t(F1) %*% Q %*% A.c - lambda*t(F1) %*% solve(W) %*% A.c)
    A = Q %*% A.c + T %*% A.d
    gcv = c(gcv, sum((obs-f0.smooth.hat)^2) / (n-sum(diag(A)))^2)
}


plot(f0.smooth.hat, type='l', bty='n', x=tt, ylim=range(f0(xx)))
par(new=TRUE)
#plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
plot(obs, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
par(new=TRUE)
plot(f0.hat, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=range(f0(xx)), type='l')
par(new=TRUE)
plot(f0(xx), x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=2, ylim=range(f0(xx)), type='l')


