empirical.bias.0 = list()
empirical.bias.1 = list()
empirical.bias.2 = list()
empirical.bias = list()

for (b in 1:B1) {
    #Empirical second derivative of beta.0
    sl = coefs[[b]][,4]
    cc.0 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.0 = rbind(cc.0, m2$coef)
    }
    
    empirical.bias.0[[b]] = cc.0[,2] * sapply(m[['fits']], function(x) x[['bw']]) / 10
    
    #Empirical second derivative of beta.1
    sl = coefs[[b]][,5]
    cc.1 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.1 = rbind(cc.1, m2$coef)
    }
    
    empirical.bias.1[[b]] = cc.1[,2] * sapply(m[['fits']], function(x) x[['bw']]) / 10
    
    #Empirical second derivative of beta.2
    sl = coefs[[b]][,6]
    cc.2 = matrix(0,0,2)
    for (i in 1:n) {
        h = m$fits[[i]]$bw
        w = epanechnikov(abs(tt-tt[i]), h)
        X = tt - tt[i]
        
        m2 = lm(sl~X, weights=w)
        cc.2 = rbind(cc.2, m2$coef)
    }
    
    empirical.bias.2[[b]] = cc.2[,2] * sapply(m[['fits']], function(x) x[['bw']]) / 10
    
    empirical.bias[[b]] = cbind(empirical.bias.0[[b]], empirical.bias.1[[b]], empirical.bias.2[[b]])
}


