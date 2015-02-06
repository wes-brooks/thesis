valid = (1:20)[c(-7,-11,-13,-16)]

list(rep(0,400))

coef = list(B1, B2)
smooth = list(B1.smooth, B2.smooth)

for (m in 1:length(test.loc)) {
    pdf(paste("~/Desktop/coverage", m, "pdf", sep="."), width=12, height=8)
    par(mar=c(2.1,2.1,1.1,1.1))
    layout(matrix(1:2,2,1))
    for (j in 1:2) {
        yy = range(sapply(valid, function(k) res.CI[[k]][[j+1]][test.loc[m],]),
                   coef[[j]][test.loc][m],
                   sapply(valid, function(k) smooth[[j]][[k]][test.loc][m]))
        
        plot.new()
        plot.window(xlim=c(0, length(valid)+1), ylim=yy)
        
        axis(2)
        
        #Put the true coefficients on the plot:
        #points(x=0, y=B1[test.loc][1])
        abline(h=coef[[j]][test.loc][m], lty=2)
        
        #Put the smoothed coefficients on the plot:
        points(x=1:length(valid)+0.1, y=sapply(valid, function(k) smooth[[j]][[k]][test.loc][m]), pch=20, cex=2, col='red')
        
        for (i in 1:length(valid)) {
            k = valid[i]
            lines(x=rep(i+0.2,2), y=res.CI[[k]][[j+1]][test.loc[m],], col=grey(0.1), lwd=2)
            lines(x=rep(i+0.35,2), y=bayes.CI[[k]][[j+1]][test.loc[m],], col=grey(0.6), lwd=2)
            lines(x=rep(i+0.5,2), y=orc.CI[[k]][[j+1]][test.loc[m],], col=grey(0.2), lwd=2)
        }
    }
    dev.off()
}