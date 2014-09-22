#SUMMARIZE RESULTS:
require(xtable)
require(brooks)
load("~/git/gwr/scratch/boston-preset-bw.RData")

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




#PLOT OF ESTIMATED COEFFICIENTS:
require(brooks)
load("~/git/gwr/scratch/boston-preset-bw.RData")

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

pdf("~/git/gwr/writeup/estimation-standalone/figure/boston-plots.pdf", 8, 8)
multiplot(plotlist=bmap, cols=2)    
dev.off()


#SUMMARY TABLE:
require(xtable)
require(brooks)
load("~/git/gwr/scratch/boston-preset-bw.RData")

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

#Generate the table, caption it, and print it with emphasis.
boston.coef.table = xtable(boston.coef.summary, align="|c|ccc|")
caption(boston.coef.table) = "The mean, standard deviation, and proportion of zeros among the estimates of the local coefficients in a model for the median house price in census tracts in Boston, with coefficients selected and fitted by local adaptive grouped regularization. The covariates are CRIM (per capita crime rate in the census tract), RM (average number of rooms per home sold in the census tract), RAD (an index of the tract's access to radial roads), TAX (property tax per USD10,000 of property value), and LSTAT (percentage of the tract's residents who are considered ``lower status\")."
label(boston.coef.table) = "tab:boston-coefs-lagr"
print(boston.coef.table, table.placement=NULL)