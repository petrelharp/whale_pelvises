#!/usr/bin/R

load("mcmc-setup.RData")
load("thedata-and-covmatrices.Rdata")

show( load("mcmcs/mcmc-run-3937.RData") )

nvars <- length(initpar)
havedata <- !is.na(thedata)
datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data

mcrun.files <- list.files('mcmcs/',"mcmc-run-.*RData", full.names=TRUE)
mcruns <- lapply( mcrun.files, function (x) { tmpenv <- environment(); tmp <- load(x,envir=tmpenv); names(tmp) <- tmp; lapply( tmp, get, envir=tmpenv ) } )

lls <- rep(NA,length(mcruns))
for (k in 1:length(mcruns)) { lls[k] <- llfun( mcruns[[k]]$mcrun$final ) }
luds <- rep(NA,length(mcruns))
for (k in 1:length(mcruns)) { luds[k] <- mcruns[[k]]$mcrun$lud( mcruns[[k]]$mcrun$final ) }

####
# traces
layout( matrix(1:(2*length(mcruns)),nrow=2,byrow=TRUE) )
lapply( mcruns, function (x) matplot(x$mcrun$batch,type='l') )
lapply( mcruns, function (x) matplot(x$mcrun$batch,type='l',log='y') )


###
# Look at residuals
fullmats <- lapply( mcruns, function (x) make.fullmat( x$mcrun$final ) )
submats <- lapply( fullmats, function (fullmat) { ((crossprod(pmat, fullmat[havedata,havedata]) %*% pmat)) } )
chols <- lapply( submats, chol )
normresids <- sapply( chols, function (fchol) { backsolve( fchol, datavec, transpose=TRUE ) } )
init.fullmat <- make.fullmat( initpar )
init.normresid <- backsolve( chol( crossprod(pmat, init.fullmat[havedata,havedata]) %*% pmat ), datavec, transpose=TRUE )

layout(1)
matplot( normresids[ order(rowSums(normresids)), ], type='l' )

####
# Looking normal yet?
qqnorms <- apply( normresids, 2, qqnorm, plot.it=FALSE )
plot( 0, type='n', xlim=range(lapply(qqnorms,"[[",'x')), ylim=range(lapply(qqnorms,"[[",'y')), xlab='', ylab='' )
lapply( seq_along(qqnorms), function (k) points(qqnorms[[k]], col=rainbow(10)[k]) )
points( qqnorm( init.normresid, plot.it=FALSE ), pch=20 )
abline(0,1)


####
# Take this one

load("mcmcs/mcmc-run-3937.RData")
load("mcmcs/mcmc-run-6418.RData"
estpar <- colMeans( mcrun$batch[1e4 + (1:(nrow(mcrun$batch)-1e4)),] )
pars <- as.data.frame( rbind(initpar,estpar) )
pars$deltaT <- pars$sigmaLdeltaT / pars$sigmaL
pars$deltaP <- pars$sigmaLdeltaP / pars$sigmaL
pars$deltaR <- pars$sigmaLdeltaR / pars$sigmaL
pars$deltaP <- pars$zetaLdeltaP / pars$zetaL
pars$deltaR <- pars$zetaLdeltaR / pars$zetaL

layout( matrix(1:nvars^2,nrow=nvars) )
opar <- par(mar=c(0,0,0,0)+.1)
for (j in 1:nvars) for (k in 1:nvars) {
    if (j==k) { 
        plot(0,type='n',xlim=c(-1,1),ylim=c(-1,1)); text(0,0,names(mcrun$initial)[j]) 
    } else {
        plot( mcrun$batch[,j], mcrun$batch[,k], pch=20, col=adjustcolor(rainbow(64),.1)[ceiling(64*(1:nrow(mcrun$batch))/(nrow(mcrun$batch+1)))] )
    }
}
par(opar)

