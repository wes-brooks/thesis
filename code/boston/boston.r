require(lagr)
require(sp)
require(doMC)
registerDoMC(cores=3)

#The data is from package 'spdep'
require(spdep)
data(boston)
boston = boston.c
boston$CHAS = as.numeric(boston$CHAS)

bw.boston = list()
bw.boston = lagr.tune(MEDV~CRIM+RM+RAD+TAX+LSTAT, data=boston, family=gaussian, coords=c('LON','LAT'), longlat=TRUE, varselect.method="wAIC", kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", verbose=TRUE)

#bw.boston[['bw']] = 0.1
#Rprof("~/Desktop/boston-profile.txt")
m.boston = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT, data=boston, family=gaussian, coords=c('LON','LAT'), fit.loc=boston[38,c("LON","LAT")], longlat=TRUE, varselect.method='wAIC', kernel=epanechnikov, bw=bw.boston[['bw']], bw.type='knn', verbose=TRUE, family='gaussian', resid.type='pearson')
#Rprof(NULL)
