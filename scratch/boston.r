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
gpclibPermit()
require(spdep)
data(boston)

bTR <- boston.c$TRACT

## MA 1970 tracts
MAtr70 <- readOGR(".", "MA-1970-tracts")

## counties in the BOSTON SMSA
MAtr70$TRACTBASE = as.integer(substring(as.character(MAtr70[['GISJOIN']]),9,12))
names(MAtr70)[3] = "ID"

## reorder data set
mm <- match(MAtr70$TRACTBASE, boston.c$TRACT)
boston.tracts = MAtr70[!is.na(mm),]
boston.map = poly_coords(boston.tracts)

#Longitude range
rlon = range(boston.map$PolyCoordsY)
# -71.52262 -70.60876

#Latitude range
rlat = range(boston.map$PolyCoordsX)
# 42.00315 42.67316


qplot(PolyCoordsY, PolyCoordsX, data=boston.map, geom="polygon", group=Poly_Name)

#All downtown Boston tracts:
dtn.tracts = MAtr70[which(MAtr70$TRACTBASE<1000),]
dtn.map = poly_coords(dtn.tracts)
qplot(PolyCoordsY, PolyCoordsX, data=dtn.map, geom="polygon", group=Poly_Name)

