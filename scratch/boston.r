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

source("scratch/interpolate.bw.r")

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

#Find the optimal bandwidth:
bw.boston = lagr.sel(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston.c, coords=boston.c[,c('LON','LAT')], longlat=TRUE, varselect.method="AICc", range=c(0,1), kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method="AICc", verbose=TRUE, family='gaussian', resid.type='pearson')
bws = interpolate.bw(bw.boston[['trace']], S=20)

#save the bandwidth object:
save(bw.boston, file="bw.boston.lagr.RData")

#Make a lagr model for the optimal bandwidth and save it:
i=0
model = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boston.c, coords=boston.c[,c('LON','LAT')], fit.loc=boston.loc, longlat=TRUE, varselect.method='AICc', kernel=epanechnikov, bw=bw.boston[['bw']], bw.type='knn', verbose=TRUE, family='gaussian', resid.type='pearson')
save(model, file=paste("boston.model.", i, ".RData", sep=""))
rm(model)
gc()

for (bw in bws) {
    i = i+1
    print(paste("iteration: ", i, ", bw: ", round(bw, 4), sep=""))
    indx = sample(1:nrow(boston.c), replace=TRUE)
    boot = boston.c[indx,]
    model = lagr(MEDV~CRIM+RM+RAD+TAX+LSTAT-1, data=boot, coords=boot[,c('LON','LAT')], fit.loc=boston.loc, longlat=TRUE, varselect.method='AICc', kernel=epanechnikov, bw=bw, bw.type='knn', verbose=TRUE, family='gaussian', resid.type='pearson')
    
    #save the bootstrapped model:
    save(model, file=paste("boston.model.", i, ".RData", sep=""))
    rm(model)
    gc()
}



#for (v in c('CRIM', 'RM', 'RAD', 'TAX', 'LSTAT')) {
#    boston.tracts@data[[paste('coef', v, sep='')]] = sapply(m.boston[['model']][['models']], function(x) x[['coef']][[v]])
#}

#Draw a map:
#boston.map = poly_coords(boston.tracts)
#qplot(PolyCoordsY, PolyCoordsX, data=boston.map, geom="polygon", group=Poly_Name, fill=coefLSTAT)
