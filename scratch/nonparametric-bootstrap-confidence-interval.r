S = 45

        cc = matrix(0, 0, 5)

for (iter in 1:370) {
    i = ((iter-1) %/% S) + 1
    j = ((iter-1) %% S) + 1

    if (j==20) {


        row=1
            a = read.table(paste(paste("output", iter, "aic", sep="/"), "jacknife", i, j, row, "txt", sep="."))
            b = read.table(paste(paste("output", iter, "beta", sep="/"), "jacknife", i, j, row, "txt", sep="."))
            cc = rbind(cc, b[,which.min(a)])
        
    }
}


t. = colMeans(cc)


tt = list()
Y.delta = matrix(0,196,0)
for (i in 1:8) {
    tt[[i]] = t(cc[i, , drop=FALSE] - t.)
    Y.delta = cbind(Y.delta, Yij[,i]-Y.j)
}

cov = list()
for (m in 1:196) {
    cov[[m]] = 0
    for (i in 1:8) {
        cov.j = Y.delta[m,i] * tt[[i]] / 8
        cov[[m]] = cov[[m]] + cov.j
    }
}

cov.2 = matrix(0, 5,5)
for (m in 1:196) {
    cov.2 = cov.2 + cov[[m]] %*% t(cov[[m]])
}