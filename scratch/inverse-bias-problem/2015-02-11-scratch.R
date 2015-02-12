#THIS ITERATIVE ALGORITHM ATTEMPTS TO MINIMIZE DIFFERENCE BOTH IN THE SMOOTHED FUNCTION AND ITS DERIVATIVE
#BUT THERE'S A PROBLEM AND THE CONDITIONAL NON-NEGATIVE DEFINITENESS IS NOT WORKING
#PROBABLY BECAUSE t(T1) %*% c ISNT ZERO.
#THIS PROBABLY MEANS I'VE MESSED UP THE MATH.

d=matrix(0,2,1)
F2 %*% solve(t(F2)%*%t(W)%*%W%*%Q0%*%F2 + diag(rep(lambda, 98))) %*% t(F2)%*%t(W)%*%(obs  - T%*%d) -> c
solve(t(F1)%*%W%*%F1%*%R) %*% t(F1) %*% t(W)%*%(obs - W%*%Q0%*%c) -> d
d


lambda=0.1
F2 %*% solve(t(F2) %*% (Q + lambda*QtI%*%Q0) %*% F2) %*% t(F2) %*% obs -> c
solve(R) %*% t(F1) %*% (obs - (Q - lambda* QtI%*%Q0) %*% c) -> d


for (i in 1:100) {
    ####Iterative: get c
    ((t(W)%*%W + t(W.)%*%W.)%*%Q0 + diag(rep(lambda,100))) %>% solve ->a2
    t(W)%*%obs + t(W.)%*%obs. - (t(W)%*%W + t(W.)%*%W.)%*%T0%*%d ->a3
    a2%*%a3 -> c.new
    #######Then d:
    solve(t(T1)%*%T1 + t(T2)%*%T2) -> tinv1
    t(T1)%*%(obs-Q1%*%c.new) + t(T2)%*%(obs.-Q2%*%c.new) ->ta
    tinv1%*% ta ->d.new
    
    #sum((d.new-d)^2 * 1) %>% print
    #sum((c.new-c)^2 * 1) %>% print
    
    c = c + (c.new-c)*1
    d = d + (d.new-d)*1
    
    print(sum((obs - (Q1 %*% c + T1 %*% d))^2 + (obs. - (Q2 %*% c + T2 %*% d))^2) + (t(c)%*%Q0%*%c)[1,1])
    print(sum((obs - (Q1 %*% c + T1 %*% d))^2 + (t(c)%*%Q0%*%c)[1,1])
}
f0.smooth.hat = (Q1 %*% c + T1 %*% d) 
f0.hat = (Q0 %*% c + T0 %*% d) 
f0..smooth.hat = Q2%*%c + T2%*%d

##############
#PROBLEM!!!
#THIS IS SUPPOSED TO BE ZERO:
#t(T1) %*% c

plot(f0.smooth.hat, type='l', bty='n', x=tt, ylim=range(f0(xx)), lwd=2)
par(new=TRUE)
#plot(f0.smooth, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l')
plot(obs, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='red', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0.hat, x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', col='blue', ylim=range(f0(xx)), type='l', lwd=2)
par(new=TRUE)
plot(f0(xx), x=tt, bty='n', ann=FALSE, xaxt='n', yaxt='n', lty=2, ylim=range(f0(xx)), type='l', lwd=2)

plot(f0..smooth.hat, type='l')
par(new=TRUE)
plot(coefs[,4], type='l', col='red', ylim=range(f0..smooth.hat))
par(new=TRUE)
plot(f0..smooth, type='l', lty=2, ylim=range(f0..smooth.hat))
