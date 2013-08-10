#!/usr/bin/R

stuff <- load("../thedata-and-covmatrices.Rdata")
true.data <- thedata
thedata <- read.csv("observed-simdata.csv",header=TRUE)
save(list=c('true.data',stuff), file="thedata-and-covmatrices.Rdata")

