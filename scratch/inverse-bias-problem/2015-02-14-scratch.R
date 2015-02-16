fn.old=Inf
while (fn.old>fn.new) {
    i=i+1
    fn.old = fn(p)
    (diag(rep(1,100)) - T0%*%solve(t(T0)%*%T0)%*%t(T0)) %*% gr(p)[1:100] -> a
    a = a / (t(a)%*%a)[1] #normalize a to length 1
    
    c = c - a/1900
    d = solve(t(T1)%*%T1 + t(T2)%*%T2) %*% (t(T1)%*%obs +t(T2)%*%obs. - t(T1)%*%Q1%*%c - t(T2)%*%Q2%*%c)
    p = c(c,d)
    fn.new=fn(p)
    cat(paste(i, ": ", fn.new, "\n", sep=""))
}