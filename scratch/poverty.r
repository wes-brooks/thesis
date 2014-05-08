#Import external libraries
require(devtools)
require(lagr)
require(SGL)
registerDoMC(cores=3)
require(RCurl)

#Import the data
source_url('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

#Establish lists to hold the bandwidths
bw = list()
bw[['lagr']] = list()

#Establish lists to hold the models
model = list()
model[['lagr']] = list()

years = c(1960, 1970, 1980, 1990, 2000, 2006)
years = c(2000)

for (yr in years) {
    #Isolate one year of data
    year = as.character(yr)
    df = pov2[pov2$year==yr,]
    
    #Define which variables we'll use as predictors of poverty:
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))

    #Lasso model
    bw[['lagr']][[year]] = lagr.sel(formula=f, data=df, family='gaussian', coords=df[,c('x','y')], longlat=TRUE, varselect.method='AICc', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', parallel=TRUE, interact=TRUE, verbose=TRUE, bwselect.method='AICc', resid.type='pearson')
    model[['lagr']][[year]] = lagr(formula=f, data=df, family='gaussian', coords=df[,c('x','y')], longlat=TRUE, N=1, varselect.method='AICc', bw=bw[['lagr']][[year]][['bw']], kernel=epanechnikov, bw.type='knn', simulation=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, resid.type='pearson')
}