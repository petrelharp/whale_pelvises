#!/usr/bin/R

load("mcmc-setup.RData")
load("thedata-and-covmatrices.Rdata")

whichsex <- rev(strsplit(getwd(),"/")[[1]])[1]
if (whichsex=="females") {
    load("mcmcs/mcmc-run-8467.RData")
} else if (whichsex=="males") {
} else { stop("wrong directory -- males or females") }

nonnegs <- c("sigmaL", "betaT", "sigmaR", "sigmaP", "zetaL", "zetaR", "omegaR", "zetaP", "omegaP" )
nonneg.inds <- match( nonnegs, names(initpar) ) 

nvars <- length(initpar)
havedata <- !is.na(thedata)
datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data

mcrun.files <- list.files('mcmcs',"mcmc-run-.*RData", full.names=TRUE)
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

if (interactive()) {
    layout(1)
    matplot( normresids[ order(rowSums(normresids)), ], type='l' )
}

####
# Looking normal yet?
if (interactive()) {
    qqnorms <- apply( normresids, 2, qqnorm, plot.it=FALSE )
    plot( 0, type='n', xlim=range(lapply(qqnorms,"[[",'x')), ylim=range(lapply(qqnorms,"[[",'y')), xlab='', ylab='' )
    lapply( seq_along(qqnorms), function (k) points(qqnorms[[k]], col=rainbow(10)[k]) )
    points( qqnorm( init.normresid, plot.it=FALSE ), pch=20 )
    abline(0,1)
}


####
# Take this one


burnin <- 2e4
usethese <- burnin + (1:(nrow(mcrun$batch)-burnin))
estpar <- colMeans( mcrun$batch[usethese,] )
pars <- as.data.frame( rbind(initpar[names(mcrun$initial)],estpar) )
quants <- as.data.frame( apply(mcrun$batch[usethese,],2,quantile,c(.05,.25,.75,.95)) )
names(quants) <- names(pars)
for ( x in c("deltaT", "deltaP","deltaR") ) {
    for (y in c("sigmaL","zetaL")) {
        combname <- paste(y,x,sep='')
        jj <- match(combname,names(mcrun$initial))
        kk <- match(y,names(mcrun$initial))
        if (!is.na(jj)) {
            outname <- paste(substr(y,1,1),x,sep='')
            quants[outname] <- quantile( mcrun$batch[usethese,jj]/mcrun$batch[usethese,kk], c(.05,.25,.75,.95) )
            pars[outname] <- pars[combname] / pars[y]
        }
    }
}
pars <- rbind( pars, quants )
pars <- pars[ ! ( grepl("delta",names(pars)) & ( grepl("sigma",names(pars)) | grepl("zeta",names(pars)) ) ) ]

samples <- mcrun$batch[ seq(burnin,nrow(mcrun$batch),length.out=1e3), ]
colnames(samples) <- names(mcrun$initial)
save( pars, samples, file="results.RData" )

####
# output tables

require(xtable)
xtable( pars[c(3,4,2,5,6), grep("delta",names(pars),invert=TRUE) ], digits=3 )
xtable( pars[c(3,4,2,5,6), grep("delta",names(pars)) ], digits=3 )

names(estpar) <- names(mcrun$initial)
species.params <- estpar[species.parindices]
sample.params <- estpar[sample.parindices]
species.transmat[species.Tind] <- species.params[species.Sind]
sample.transmat[sample.Tind] <- sample.params[sample.Sind] * sample.Pcoef
species.covmat <- as.matrix( tcrossprod(species.transmat) )
sample.covmat <- as.matrix( tcrossprod(sample.transmat) )
species.subcovmat <- as.matrix( tcrossprod(species.transmat[-1,-1]) )
sample.subcovmat <- as.matrix( tcrossprod(sample.transmat[-1,-1]) )
species.notestes.covmat <- as.matrix( tcrossprod(species.transmat[-2,-2]) )
sample.notestes.covmat <- as.matrix( tcrossprod(sample.transmat[-2,-2]) )

# how much faster do pelvic bones evolve with testes doing their thing?
species.covmat[5,5] / species.notestes.covmat[4,4]
# both: 1.280029
# females: 1.291685

# correlation matrices
xtable( cov2cor(species.covmat)[c(1,2,3,5),c(1,2,3,5)], digits=2 )
xtable( cov2cor(sample.covmat)[c(1,3:6),c(1,3:6)], digits=2 )
xtable( cov2cor(species.subcovmat[c(1,2,4),c(1,2,4)]) )

# get posterior distribution on the correlations
get.correlations <- function (par) {
    species.params <- par[species.parindices]
    sample.params <- par[sample.parindices]
    species.transmat[species.Tind] <- species.params[species.Sind]
    sample.transmat[sample.Tind] <- sample.params[sample.Sind] * sample.Pcoef
    species.subcovmat <- as.matrix( tcrossprod(species.transmat[-1,-1]) )
    sample.subcovmat <- as.matrix( tcrossprod(sample.transmat[-1,-1]) )
    return( cov2cor(species.subcovmat[c(1,2,4),c(1,2,4)]) )
}
posterior.cors <- apply( samples, 1, get.correlations )
dim( posterior.cors ) <- c(3,3,nrow(samples))

pdf(file="posterior-correlations.pdf", width=4, height=3, pointsize=10 )
par(mar=c(4,4,1,1)+.1)
hist( posterior.cors[1,2,], breaks=30, xlim=range(posterior.cors), ylim=c(0, 6), col=adjustcolor('black',.5), xlab="value", main='', ylab="posterior density", freq=FALSE )
hist( posterior.cors[1,3,], breaks=30, col=adjustcolor('red',.5), add=TRUE, freq=FALSE )
hist( posterior.cors[2,3,], breaks=30, col=adjustcolor('blue',.5), add=TRUE, freq=FALSE )
legend("topright", fill=adjustcolor(c('black','red','blue'),.5), legend=c("testes-ribs", "testes-pelvis", "ribs-pelvis"), title="correlations" )
dev.off()

cors.df <-  data.frame( 'testes-ribs'=posterior.cors[1,2,], 'testes-pelvis'=posterior.cors[1,3,], 'ribs-pelvis'=posterior.cors[2,3,] )
rbind( sapply( cors.df, quantile, prob=c(.025,.975)), sapply( cors.df, summary ) )
#         testes.ribs testes.pelvis ribs.pelvis
# 2.5%     -0.5327032     0.1788103  -0.3354422
# 97.5%     0.3969590     0.8008125   0.2472017
# Min.       -0.69880       -0.1143    -0.58270
# 1st Qu.    -0.28530        0.4401    -0.14800
# Median     -0.09575        0.5697    -0.04419
# Mean       -0.09502        0.5457    -0.05181
# 3rd Qu.     0.07716        0.6741     0.03618
# Max.        0.63000        0.8752     0.48420

###
# pairwise correlations
if (interactive()) {
    layout( matrix(1:nvars^2,nrow=nvars) )
    opar <- par(mar=c(0,0,0,0)+.1)
    subsamp <- floor( seq( 1, nrow(mcrun$batch), length.out=1000 ) )
    for (j in 1:nvars) for (k in 1:nvars) {
        if (j==k) { 
            plot(0,type='n',xlim=c(-1,1),ylim=c(-1,1)); text(0,0,names(mcrun$initial)[j]) 
        } else {
            plot( mcrun$batch[subsamp,j], mcrun$batch[subsamp,k], pch=20, col=adjustcolor(rainbow(64),.1)[ceiling(64*subsamp/(nrow(mcrun$batch+1)))] )
        }
    }
    par(opar)
}

