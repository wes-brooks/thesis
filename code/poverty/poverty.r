#Import external libraries
require(devtools)
install_github('SGL', 'wrbrooks')
install_github('gwselect', 'wesesque', ref='agLasso')
require(SGL)
require(gwselect)
#require(spgwr)
registerDoMC(cores=3)
require(RCurl)

#Import the data
source_url('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

#Establish lists to hold the bandwidths
bw = list()
bw[['GWAL']] = list()
bw[['GWR']] = list()

#Establish lists to hold the models
model = list()
model[['GWAL']] = list()
model[['GWR']] = list()

years = c(1960, 1970, 1980, 1990, 2000, 2006)
years = c(1970)

for (yr in years) {
    #Isolate one year of data
    year = as.character(yr)
    df = pov2[pov2$year==yr,]
    
    ####Produce the models via lasso and via elastic net:
    #Define which variables we'll use as predictors of poverty:
    #predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
    f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))

    #Lasso model
    p = lagr(logitindpov~pag, data=df, coords=c('x','y'), varselect.method='wAIC', kernel=epanechnikov, bw.type='knn', bw=0.2, verbose=TRUE)
    bw = lagr.tune(formula=f, data=df, family='gaussian', coords=c('x','y'), varselect.method='wAIC', kernel=epanechnikov, bw.type='knn', bwselect.method='AICc', verbose=TRUE, n.lambda=80, lagr.convergence.tol=0.005, lambda.min.ratio=0.01)    
    bw[['GWAL']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', coords=df[,c('x','y')], longlat=TRUE, varselect.method='AICc', gweight=epanechnikov, tol.bw=0.01, bw.type='knn', parallel=TRUE, interact=TRUE, verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
    model[['GWAL']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', coords=df[,c('x','y')], longlat=TRUE, N=1, varselect.method='AICc', bw=bw[['GWAL']][[year]][['bw']], gweight=epanechnikov, bw.type='knn', simulation=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, resid.type='pearson')

    #Use my code to do the traditional GWR; (currently buggy)
    oracle = lapply(1:533, function(x) {return(predictors)})
    bw[['GWR']][[year]] = gwglmnet.sel(f, data=df, oracle=oracle, coords=df[,c('x','y')], family='gaussian', mode.select='AICc', longlat=TRUE, gweight=epanechnikov, tol.bw=0.01, bw.method='knn', parallel=TRUE, interact=FALSE, verbose=TRUE, bw.select='AICc', resid.type='pearson')
    model[['GWR']][[year]] = gwglmnet(f, data=df, oracle=oracle, coords=df[,c('x','y')], family='gaussian', mode.select='AICc', longlat=TRUE, bw=bw[['GWR']][[year]][['bw']], gweight=epanechnikov, bw.method='knn', simulation=TRUE, parallel=TRUE, interact=FALSE, verbose=TRUE)
}