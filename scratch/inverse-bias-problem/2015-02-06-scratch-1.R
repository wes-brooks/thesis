inv = vector()
for (i in 1:400) {
    h = model$fits[[i]]$bw
    dist = D[i,]
    wt = lagr:::epanechnikov(dist, h)
    X.loc = cbind(1, location$x-location$x[i], location$y-location$y[i])
    W = diag(wt)
    
    P = ((X.loc %*% solve(t(X.loc) %*% W %*% X.loc) %*% t(X.loc) %*% W))
    Pp = P[which(P[i,]!=0),]
    inv = c(inv, (ginv(t(Pp) %*% Pp) %*% t(P) %*% B1.smooth[[20]])[i,1])
}