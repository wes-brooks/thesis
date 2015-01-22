library(lagr)
## Conversion table 1980/1970
# ICPSR_07913.zip
# 07913-0001-Data.txt
# http://dx.doi.org/10.3886/ICPSR07913.v1
# Provider: ICPSR
# Content: text/plain; charset="us-ascii"
#
# TY - DATA
# T1 - Census of Population and Housing 1980 [United States]:
# 1970-Pre 1980 Tract Relationships
# AU - United States Department of Commerce. Bureau of the Censusboston 21
# DO - 10.3886/ICPSR07913.v1
# PY - 1984-06-28
# UR - http://dx.doi.org/10.3886/ICPSR07913.v1
# PB - Inter-university Consortium for Political and Social Research
# (ICPSR) [distributor]
# ER -

poly_coords<- function(shapefile) {
    if (nrow(data.frame(shapefile$ID)) < 1) {
        print ("No ID field in SpatialPolygon")
    } else {
        Order <-0 
        YX3 <- as.numeric("XX", "XX", "XX", "XX")
        num_polys <- nrow(shapefile@data)+1
        YX3 <- as.numeric("XX", "XX", "XX")
        
        curr_poly <- shapefile@data[1,]
        curr_poly_start_row <- 1
        
        
        for (curr_row in curr_poly_start_row:num_polys) {
            curr_poly_row<-shapefile@data[curr_row,]
            curr_poly_end_row = curr_row - 1	
            Poly_n= shapefile@data[curr_poly_start_row:curr_poly_end_row,]
            curr_poly_start_row = curr_row
            Poly_Name<-as.vector(Poly_n$ID)
            Poly<-shapefile[shapefile$ID==Poly_Name,]
            PolyCoords<-lapply(slot(Poly, "polygons"), function(x) lapply(slot(x,"Polygons"), function(y) slot(y,"coords")))
            PolyCoordsY<-PolyCoords[[1]][[1]][,1]
            PolyCoordsX<-PolyCoords[[1]][[1]][,2]
            Order<- 1:nrow(data.frame(PolyCoordsX)) + max(Order)
            YX1<- data.frame(Poly_Name, Order, PolyCoordsY, PolyCoordsX)
            YX2<-rbind(YX3,YX1)
            YX3<-YX2
        }
        join<-merge(YX3, shapefile@data, by.x="Poly_Name", by.y= "ID", all=T)
        join[order(join$Order),][1:nrow(join)-1,]
    }
}

## match against boston data set
require(ggplot2)
require(rgeos)
require(rgdal)
require(maptools)
require(spdep)
require(lagr)
require(doMC)
registerDoMC(7)

require(xtable)
require(brooks)

#source("scratch/interpolate.bw.r")

#Import the boston house price dataset
data(boston)
boston.c$CHAS = as.numeric(boston.c$CHAS)

## MA 1970 tracts
setwd("data/1970-tract-shapes-MA/")
MAtr70 <- readOGR(".", "MA-1970-tracts")
setwd("../..")

## counties in the BOSTON SMSA
MAtr70$TRACTBASE = as.integer(substring(as.character(MAtr70[['GISJOIN']]),9,12))
names(MAtr70)[3] = "ID"

## reorder data set
mm <- match(MAtr70$TRACTBASE, boston.c$TRACT)
boston.tracts = MAtr70[!is.na(mm),]

#All downtown Boston tracts:
dtn.tracts = MAtr70[which(MAtr70$TRACTBASE<1000),]

#Add the missing tracts:
indx = which(!dtn.tracts$TRACTBASE %in% boston.tracts$TRACTBASE)
boston.tracts = rbind(boston.tracts, dtn.tracts[indx,])

#
dtn.loc = as.data.frame(t(sapply(indx, function(i) dtn.tracts@polygons[[i]]@labpt)))
colnames(dtn.loc) = c("LON", "LAT")

#Create a list of locations
boston.loc = rbind(boston.c[,c("LON", "LAT")], dtn.loc)
rownames(boston.loc) = boston.tracts@data$TRACTBASE

#Make a lagr model for the bandwidth and save it:
bw = lagr.tune(MEDV~CRIM+RM+RAD+TAX+LSTAT, data=bd, coords=c('LON','LAT'), verbose=TRUE, longlat=TRUE, varselect.method='wAICc', bwselect.method='AICc', tol.loc=0.01, kernel=epanechnikov, bw.type='knn', lagr.convergence.tol=0.05, n.lambda=200, lambda.min.ratio=1e-7, lagr.max.iter=50, family='gaussian')
model = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT, data=boston.c, coords=c('LON','LAT'), fit.loc=boston.loc, longlat=TRUE, varselect.method='AIC', kernel=epanechnikov, bw=0.26, bw.type='knn', verbose=TRUE, lagr.convergence.tol=0.05, n.lambda=200, lambda.min.ratio=1e-7, lagr.max.iter=50, family='gaussian')

cc = as.data.frame(t(sapply(model[['model']], function(x) x[['coef']])))


vars = c('CRIM', 'RM', 'RAD', 'TAX', 'LSTAT')

boston.coef.summary = matrix(NA, nrow=0, ncol=3)
for (v in vars) {
    row = vector()
    colname.coef = paste("coef", v, sep="")
    
    #Add the table's elements to this row:
    row = c(row, mean(boston.tracts@data[,colname.coef]))
    row = c(row, sd(boston.tracts@data[,colname.coef]))
    row = c(row, mean(boston.tracts@data[,colname.coef]==0))
    
    #Add this row to the table:
    boston.coef.summary = rbind(boston.coef.summary, row)
}
rownames(boston.coef.summary) = vars
colnames(boston.coef.summary) = c('Mean', 'SD', 'Prop. zero')

bb = as.factor(apply(sapply(model[['model']], function(x) ifelse(x[['coef']][-1]==0,0,1)), 2, function(x) paste(x,collapse="")))


#PLOT OF ESTIMATED COEFFICIENTS:
for (v in c('CRIM', 'RM', 'RAD', 'TAX', 'LSTAT', '(Intercept)')) {
    boston.tracts@data[[paste('coef', v, sep='')]] = sapply(model[['model']], function(x) x[['coef']][[v]])
}
boston.tracts@data[['model']] = bb

#Draw a map:
boston.map = poly_coords(boston.tracts)
bmap = list()
for (v in c('CRIM', 'RM', 'RAD', 'TAX', 'LSTAT')) {
    bmap[[v]] = ggplot(boston.map) +
        aes(x=PolyCoordsY, y=PolyCoordsX, group=Poly_Name) +
        aes_string(fill=paste('coef', v, sep='')) +
        geom_polygon() +
        scale_fill_gradient2(low='orange', mid='white', high="purple", midpoint=0) +
        xlab("longitude") +
        ylab("latitude")
}
bmap[['model']] = ggplot(boston.map) +
        aes(x=PolyCoordsY, y=PolyCoordsX, group=Poly_Name) +
        aes(fill=model) +
        geom_polygon() +
        scale_fill_brewer(type='qual', labels=c("RM, LSTAT", "RM, RAD, LSTAT", "CRIM, RM,\nRAD, LSTAT", "CRIM, RM,\nRAD,TAX, LSTAT")) +
        theme(legend.text=element_text(size=rel(0.6))) +
        xlab("longitude") +
        ylab("latitude")

pdf("~/git/gwr/writeup/estimation-standalone/figure/boston-plots.pdf", 7, 8)
multiplot(plotlist=bmap, cols=2)    
dev.off()


#SUMMARY TABLE:
boston.coef.summary = matrix(NA, nrow=0, ncol=3)
for (v in c('CRIM', 'RM', 'RAD', 'TAX', 'LSTAT')) {
    row = vector()
    colname.coef = paste("coef", v, sep="")
    
    #Add the table's elements to this row:
    row = c(row, mean(boston.tracts@data[1:506,colname.coef]))
    row = c(row, sd(boston.tracts@data[1:506,colname.coef]))
    row = c(row, sum(boston.tracts@data[1:506,colname.coef]==0))
    
    #Add this row to the table:
    boston.coef.summary = rbind(boston.coef.summary, row)
}
rownames(boston.coef.summary) = vars
colnames(boston.coef.summary) = c('Mean', 'SD', '\\begin{tabular}{c}Zero coef. \\\\ count\\end{tabular}')

#Generate the table, caption it, and print it with emphasis.
boston.coef.table = xtable(boston.coef.summary, align="|c|ccc|", digits=c(2,2,2,0))
caption(boston.coef.table) = "The mean, standard deviation, and proportion of zeros among the estimates of the local coefficients in a model for the median house price in census tracts in Boston, with coefficients selected and fitted by local adaptive grouped regularization. The covariates are CRIM (per capita crime rate in the census tract), RM (average number of rooms per home sold in the census tract), RAD (an index of the tract's access to radial roads), TAX (property tax per USD10,000 of property value), and LSTAT (percentage of the tract's residents who are considered ``lower status\")."
label(boston.coef.table) = "tab:boston-coefs-lagr"
print(boston.coef.table, table.placement=NULL)

