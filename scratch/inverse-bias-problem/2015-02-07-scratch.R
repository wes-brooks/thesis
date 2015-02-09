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
    if (d %% 2 == 0) {
        res = theta(m,d) * r^(2*m-d) * log(r)
    }
    
    #odd:
    if (d %% 2 ==1) {
        res = theta(m,d) * r^(2*m-d)
    }
    
    res
}

#Thin plate smoothing spline penalty matrix:
K = matrix(0, 100, 100)
for (i in 1:100) {
    K[i,] = eta(abs(tt-tt[i]), 2, 1)
}

eta.derivative = function(t, m, d) {
#     #d=2:
#     if (d %% 2 == 0) {
#         res = r^(2*m-d) * log(r)
#     }
    
    #d = 1:
    if (d ==1) {
        res = sign(t) * (2*m-d) * abs(t)^(2*m-d-1)
    }
}

L.1 = function(t0, s, m=2, d=1, wt) {
    #The response for the inverse-bias problem is the locally linear estimate, smoothed over the kernel window
    wt*eta(0, m, d) + wt*eta.derivative(s-t0, m, d)
}

W = matrix(0, n, n)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    W[i,] = w/sum(w)
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


E. = matrix(0, n, n)
for (i in 1:n) {
    h = m$fits[[i]]$bw
    w = epanechnikov(abs(tt-tt[i]), h)
    sumw = sum(w)
    for (j in 1:n) {
        if (w[j] > 0)
            E.[i,] = E.[i,] + w[j]/sumw*eta.derivative(abs(tt-tt[j]), m=2, d=1)
    }
}

xx = seq(0,5,len=100)
F1 = qr(T) %>% qr.Q(complete=TRUE)[,1:2]
F2 = qr(T) %>% qr.Q(complete=TRUE)[,3:100]
R = qr(T) %>% qr.R

lambda=0.01
#F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% ginv(Q) %*% Q0 %*% F2) %*% t(F2) %*% f0.smooth -> c
#solve(R) %*% (t(F1) %*% f0.smooth - t(F1) %*% Q %*% c - lambda*t(F1) %*% ginv(Q) %*% Q0 %*%c) -> d
F2 %*% solve(t(F2) %*% Q %*% F2 + lambda * t(F2) %*% solve(W) %*% F2) %*% t(F2) %*% f0.smooth -> c
solve(R) %*% (t(F1) %*% f0.smooth - t(F1) %*% Q %*% c - lambda*t(F1) %*% solve(W) %*%c) -> d


f0.smooth.hat = (Q %*% c + T %*% d) 
f0.hat = (Q0 %*% c + T0 %*% d) 
plot(f0.smooth.hat, type='l', bty='n', x=tt, ylim=range(f0(xx)))
par(new=TRUE)
plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
par(new=TRUE)
plot(f0.hat, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=range(f0(xx)), type='l')
par(new=TRUE)
plot(f0(xx), x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=2, ylim=range(f0(xx)), type='l')

