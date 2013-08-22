#!/usr/bin/R

load("mcmc-setup.RData")
load("thedata-and-covmatrices.Rdata")

show( load("mcmc-run-9430.RData") )

nvars <- length(initpar)
havedata <- !is.na(thedata)
datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data

mcrun.files <- list.files('mcmcs/',"mcmc-run-.*RData")
mcruns <- lapply( mcrun.files, function (x) { tmpenv <- environment(); tmp <- load(x,envir=tmpenv); names(tmp) <- tmp; lapply( tmp, get, envir=tmpenv ) } )

lls <- rep(NA,length(mcruns))
for (k in 1:length(mcruns)) { lls[k] <- llfun( mcruns[[k]]$mcrun$final ) }
luds <- rep(NA,length(mcruns))
for (k in 1:length(mcruns)) { luds[k] <- mcruns[[k]]$mcrun$lud( mcruns[[k]]$mcrun$final ) }

####
# traces
layout( t(seq_along(mcruns)) )
lapply( mcruns, function (x) matplot(x$mcrun$batch,type='l') )


###
# Look at residuals
fullmats <- lapply( mcruns, function (x) make.fullmat( x$mcrun$final ) )
submats <- lapply( fullmats, function (fullmat) { ((crossprod(pmat, fullmat[havedata,havedata]) %*% pmat)) } )
chols <- lapply( submats, chol )
normresids <- sapply( chols, function (fchol) { backsolve( fchol, datavec, transpose=TRUE ) } )
init.fullmat <- make.fullmat( initpar )
init.normresid <- backsolve( chol( crossprod(pmat, init.fullmat[havedata,havedata]) %*% pmat ), datavec, transpose=TRUE )


matplot( normresids[ order(rowSums(normresids)), ], type='l' )

####
# Looking normal yet?
qqnorms <- apply( normresids, 2, qqnorm, plot.it=FALSE )
plot( 0, type='n', xlim=range(lapply(qqnorms,"[[",'x')), ylim=range(lapply(qqnorms,"[[",'y')), xlab='', ylab='' )
lapply( seq_along(qqnorms), function (k) points(qqnorms[[k]], col=rainbow(10)[k]) )
points( qqnorm( init.normresid, plot.it=FALSE ), pch=20 )
abline(0,1)
