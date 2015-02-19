i=0
fn.old=Inf
norma = Inf
while (norma > 0.01) {
    i=i+1
    fn.old = fn(p)
    (diag(rep(1,400)) - T0%*%solve(t(T0)%*%T0)%*%t(T0)) %*% gr(p)[1:400] -> a
    norma = (t(a)%*%a)[1]
    a = a / norma #normalize a to length 1
    a*dir.der2(a,p)[1]^(-1) -> aa
    c - aa -> c
    
    d = solve(t(T1)%*%T1 + t(T2)%*%T2 + t(T3)%*%T3) %*% (t(T1)%*%obs + t(T2)%*%obs.1 + t(T3)%*%obs.2 - t(T1)%*%Q1%*%c - t(T2)%*%Q2%*%c - t(T3)%*%Q3%*%c)
    p = c(c,d)
    fn.new=fn(p)
    cat(paste(i, ": ", fn.new, "\n", sep=""))
}

