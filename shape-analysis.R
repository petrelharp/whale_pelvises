#/usr/bin/R --vanilla

####
# 


####
## read in data
load("thedata-and-covmatrices.Rdata")
shapediff <- read.table("61_pairwise_shapes.out", header=TRUE)
