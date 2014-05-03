require(lagr)
require(sp)
require(doMC)
registerDoMC(cores=2)

load('data/boston/boston.rda')
boston = boston.c
boston$CHAS = as.numeric(boston$CHAS)

bw = list()
#bw = lagr.sel(MEDV~CRIM+INDUS+CHAS+NOX+RM+AGE+DIS+RAD+TAX+PTRATIO+B+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], longlat=TRUE, varselect.method="AICc", range=c(0,1), kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", parallel=TRUE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson')

bw[['bw']] = 0.3
Rprof("~/Desktop/boston-profile.txt")
ml = lagr(MEDV~CRIM+INDUS+CHAS+NOX+RM+AGE+DIS+RAD+TAX+PTRATIO+B+LSTAT-1, data=boston, coords=boston[,c('LON','LAT')], longlat=TRUE, N=1, varselect.method='AICc', kernel=epanechnikov, bw=bw[['bw']], bw.type='knn', parallel=FALSE, interact=TRUE, verbose=TRUE, family='gaussian', resid.type='pearson', simulation=TRUE)
Rprof(NULL)
