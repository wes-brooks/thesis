require(lagr)
require(sp)
require(doMC)
registerDoMC(cores=1)

#The data is from package 'spdep'
require(spdep)
data(boston)
boston = boston.c
boston$CHAS = as.numeric(boston$CHAS)

bw.boston = list()
bw.boston = lagr.sel(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], longlat=TRUE, varselect.method="AICc", range=c(0,.15), kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", verbose=TRUE, family='gaussian', resid.type='pearson')

bw.boston[['bw']] = 0.05
#Rprof("~/Desktop/boston-profile.txt")
m.boston = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], fit.loc=boston[1,c("LON","LAT")], longlat=TRUE, varselect.method='AICc', kernel=epanechnikov, bw=bw.boston[['bw']], bw.type='knn', verbose=TRUE, family='gaussian', resid.type='pearson')
#Rprof(NULL)
