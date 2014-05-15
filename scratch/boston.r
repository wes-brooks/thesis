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


widths <- c(ID=5L, FIPS70State=2L, FIPS70cty=3L, Tract70=6L, FIPS80State=2L, FIPS80cty=3L, f1=7L, CTC=6L, f2=2L, intersect1=3L, intersect2=3L, name=30L)
dta0 <- read.fwf("ICPSR_07913/DS0001/07913-0001-Data.txt", unname(widths), col.names=names(widths), colClasses=rep("character", 12), as.is=TRUE)
sub <- grep("25", dta0$FIPS80State)
MA <- dta0[sub,]

## match against boston data set
require(ggplot2)
require(rgeos)
require(rgdal)
require(maptools)
gpclibPermit()
require(spdep)
data(boston)

bTR <- boston.c$TRACT
x1 <- match(as.integer(MA$Tract70), bTR)
BOSTON <- MA[!is.na(x1),]

## MA 1990 tracts

MAtr90 <- readOGR(".", "MA-tracts")

## counties in the BOSTON SMSA
#BOSTON_SMSA <- MAtr90[MAtr90$CO,]
BOSTON_SMSA <- MAtr90
BOSTON_SMSA$TRACTBASE = substring(as.character(BOSTON_SMSA[['GISJOIN']]),9,12)
BOSTON_SMSA$TRACT = substring(as.character(BOSTON_SMSA[['GISJOIN']]),9)
proj4string(BOSTON_SMSA) <- CRS(paste("+proj=longlat +datum=NAD27 +no_defs", "+ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"))
CTC4 <- substring(BOSTON$CTC, 1, 4)
CTC4u <- unique(CTC4)
TB_CTC4u <- match(BOSTON_SMSA$TRACTBASE, CTC4u)

## match 1980 tracts with 1990
BOSTON_SMSA1 <- BOSTON_SMSA[!is.na(TB_CTC4u),]
## union Polygons objects with same 1970 tract code
#library(rgeos)
BOSTON_SMSA2 <- unionSpatialPolygons(BOSTON_SMSA1, IDs=as.character(BOSTON_SMSA1$TRACTBASE))
#BOSTON_SMSA2 <- gUnaryUnion(BOSTON_SMSA1, id=as.character(BOSTON_SMSA1$TRACTBASE))

## reorder data set
mm <- match(as.integer(as.character(row.names(BOSTON_SMSA2))), boston.c$TRACT)
df <- boston.c[mm,]
row.names(df) <- df$TRACT
row.names(BOSTON_SMSA2) <- as.character(as.integer(row.names(BOSTON_SMSA2)))

## create SpatialPolygonsDataFrame
BOSTON_SMSA3 <- SpatialPolygonsDataFrame(BOSTON_SMSA2, data=data.frame(poltract=row.names(BOSTON_SMSA2), row.names=row.names(BOSTON_SMSA2)))
BOSTON_SMSA4 <- spCbind(BOSTON_SMSA3, df)
mm1 <- match(boston.c$TRACT, row.names(BOSTON_SMSA4))
BOSTON_SMSA5 <- BOSTON_SMSA4[mm1,]

boston.map = poly_coords(BOSTON_SMSA5)
boston.map = boston.map[boston.map$PolyCoordsY< -70.5,]

qplot(PolyCoordsY, PolyCoordsX, data=boston.map, geom="polygon", group=Poly_Name, fill=CMEDV)

