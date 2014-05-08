require(lagr)
require(sp)
require(doMC)
registerDoMC(cores=3)

load('data/boston/boston.rda')
boston = boston.c
boston$CHAS = as.numeric(boston$CHAS)

bw.boston = list()
bw.boston = lagr.sel(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], longlat=TRUE, varselect.method="AICc", range=c(0,1), kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", parallel=TRUE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson')

#bw.boston[['bw']] = 0.2
#Rprof("~/Desktop/boston-profile.txt")
m.boston = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], longlat=TRUE, N=1, varselect.method='AICc', kernel=epanechnikov, bw=bw.boston[['bw']], bw.type='knn', parallel=FALSE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson', simulation=TRUE)
#Rprof(NULL)
