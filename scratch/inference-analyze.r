library(dplyr)
library(reshape2)

b = read.table("~/git/gwr/scratch/trace.txt", header=FALSE)
colnames(b) = c("bw", "loss")


