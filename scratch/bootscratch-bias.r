b = 1

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
