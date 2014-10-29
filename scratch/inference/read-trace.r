library(dplyr)
library(reshape2)
jack = read.table("~/Desktop/o/trace.jacknife.txt")
colnames(jack) = c("bw", "loss")

jj = acast(unique(jack), bw~., value.var='loss', fun.aggregate=mean)
