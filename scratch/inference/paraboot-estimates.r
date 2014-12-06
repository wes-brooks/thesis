# Plot bootstrap estimates:

cc2 = list()
for (i in 1:B) {
    cc2[[i]] = sapply(paraboot.models.2[[i]][['model']], function(x) x[['model']][['results']][['big.avg']])
}
cc = list()
for (i in 1:B) {
    cc[[i]] = sapply(paraboot.models.3[[i]][['model']], function(x) x[['model']][['results']][['big.avg']])
}


plot(cc[[1]][1,], type='l', ylim=c(-1,12), lty=2)
for(i in 2:B) {
    par(new=TRUE)
    plot(cc[[i]][1,], type='l', ylim=c(-1,12), lty=2)
}
par(new=TRUE)
plot(f0(tt), type='l', col='red', lwd=2, ylim=c(-1,12))

plot(cc2[[1]][1,], type='l', ylim=c(-1,12), lty=2)
for(i in 2:B) {
    par(new=TRUE)
    plot(cc2[[i]][1,], type='l', ylim=c(-1,12), lty=2)
}
par(new=TRUE)
plot(f0(tt), type='l', col='red', lwd=2, ylim=c(-1,12))





yy = c(-1,1)
l = 3

plot(cc[[1]][l,], type='l', ylim=yy, lty=2)
for(i in 2:B) {
    par(new=TRUE)
    plot(cc[[i]][l,], type='l', ylim=yy, lty=2)
}
par(new=TRUE)
plot(f2(tt), type='l', col='red', lwd=2, ylim=yy)

plot(cc2[[1]][l,], type='l', ylim=yy, lty=2)
for(i in 2:B) {
    par(new=TRUE)
    plot(cc2[[i]][l,], type='l', ylim=yy, lty=2)
}
par(new=TRUE)
plot(f2(tt), type='l', col='red', lwd=2, ylim=yy)