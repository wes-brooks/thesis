lambda=0.01

D.1 = matrix(0, 300, 304)
D.1[1:100,1:100] = Q1
D.1[1:100,101:102] = T1
D.1[101:200,103:202] = Q2
D.1[101:200,203:204] = T2
D = t(D.1) %*% D.1
D[201:300,205:304] = lambda*Q0


D.1 = matrix(0, 300, 304)
D.1[1:100,1:100] = Q0
D.1[1:100,101:102] = T0
D.1[101:200,103:202] = Q0
D.1[101:200,203:204] = T0
D = t(D.1) %*% D.1
D[201:300,205:304] = sqrt(lambda)*sqrtm(Q0)



d = cbind( cbind(t(obs), t(obs.)) %*% D.1[1:200,1:204], matrix(rep(0,100),1,100) )


#Constraint matrix:
A = matrix(0,304,204)

#Constrain all appearances of c to be equal:
a = matrix(0,304,1)
a[1] = 1
a[103] = -1
for (k in 1:100) {
    A[,k] = c(rep(0,k-1), a[1:(length(a)-k+1)])
}
for (k in 103:202) {
    A[,k-2] = c(rep(0,k-1), a[1:(length(a)-k+1)])
}

#Constrain all appearances of d to be equal:
b = matrix(0,304,1)
b[101] = 1
b[203] = -1
for (k in 1:2) {
    A[,200+k] = c(rep(0,k-1), b[1:(length(a)-k+1)])
}

#Constrain t(F1) %*% c = 0:
A[1:100,203:204] = F1


solve.QP(D, d, A, meq=204)




