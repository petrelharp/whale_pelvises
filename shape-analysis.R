#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(Matrix)
require(mcmc)

load("shape-stuff.RData")

#####
# Likelihood, species tree

# Estimation:
initpar <- c( kS=20, sigma2S=0.11, gammaP=0.02, xi2P=0.03 )
stopifnot( class(shared.paths) == "dsyMatrix"  & shared.paths@uplo == "U" )  # if so, changing upper tri also changes lower tri
make.spmat <- function ( par, edge.weights=sp.edge.testes ) {
    # par = sigma2S, gammaP
    edgelens <- ( sptree$edge.length * par[1] * exp( par[2]/par[1] * (edge.weights) ) )
    shared.paths@x[ sput ][spmap.nonz] <- as.vector( sp.mapping %*% edgelens )  # note: update @x rather than entries to preserve symmetry
    return( shared.paths )
}
thesepars <- match( c("kS", "sigma2S", "gammaP"), names(initpar) )
lud <- function (par, ...) {
    # par = ks, sigma2S, gammaP
    if (any(par[1:2]<=0)) { return( -Inf ) }
    spmat <- make.spmat( par[2:3], ... )[havedata,havedata]
    x <- ( datavec - par[1] * diag(spmat) )
    fchol <- chol(2*par[1]*spmat^2)
    return( (-1/2) * sum( (par / prior.means) ) - sum( backsolve( fchol, x, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}
# priors: exponential
prior.means <- c(20,.2,.2)
stopifnot( length(prior.means) == length(thesepars) )

have.pelvic <- !is.na(pelvic.speciesdiff[ut])
have.rib <- !is.na(rib.speciesdiff[ut])
have.both <- have.pelvic & have.rib
havedata <- have.both

###
# fix kS and look at 2D density
#  across parameter grid and sampled values in sampled.edge.testes
nsamps <- 100

ks.vals <- seq(15,23,length.out=10)
sigma2S.vals <- seq(1,5,length.out=10)
gammaP.vals <- seq(-1,2,length.out=10)
pargrid <- expand.grid( ks=ks.vals, sigma2S=sigma2S.vals, gammaP=gammaP.vals )

havedata <- have.pelvic
datavec <- pelvic.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
pargrid$pelvic <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

havedata <- have.ribs
datavec <- rib.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
pargrid$ribs <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

havedata <- ( have.pelvic & have.ribs )
datavec <- pelvic.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
pargrid$sub.pelvic <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

maxvals <- apply( pargrid[,(1+length(thesepars)):ncol(pargrid)], 2, function (x) {
        tmp <- which.max(x)
        t( as.vector( pargrid[tmp,1:length(thesepars)] ) )
    })
rownames(maxvals) <- colnames(pargrid)[1:length(thesepars)]

# plot these:
fcont <- function (z,fixed,data=pargrid,...) {
    if (is.character(z) & length(z)==1) { z <- data[[z]] }
    if (!missing(fixed)) { 
        xynames <- setdiff( names(pargrid)[1:3], fixed )
        z <- tapply( z, pargrid[xynames], max, na.rm=TRUE )
    } else {
        xynames <- c("sigma2S","gammaP")
    }
    xvals <- get(paste(xynames[1],".vals",sep=''))
    yvals <- get(paste(xynames[2],".vals",sep=''))
    dim(z) <- c(length(xvals),length(yvals))
    contour( xvals, yvals, z, levels=seq( quantile(z,.9), max(z)+1, length.out=20 ), ... )
}

fcont("pelvic", fixed="ks", col='red',xlab="sigma2S", ylab="gammaP")
fcont("ribs", fixed="ks", add=TRUE,col='blue')
fcont("sub.pelvic", fixed="ks", col='green',add=TRUE)

# effect of sampling from posterior on testes size?
require(parallel)
nreps <- 4
havedata <- have.pelvic
datavec <- pelvic.speciesdiff[ut][havedata]
pelvic.sgrids <- mclapply( 1:nreps, function (tmp) {
        k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    }, mc.cores=4 )

havedata <- have.ribs
datavec <- rib.speciesdiff[ut][havedata]
rib.sgrids <- mclapply( 1:nreps, function (tmp) {
        k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    }, mc.cores=4 )

havedata <- ( have.pelvic & have.ribs )
datavec <- pelvic.speciesdiff[ut][havedata]
sub.pelvic.sgrids <- mclapply( 1:nreps, function (tmp) {
        k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        sub.pelvic.sgrid <- apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    }, mc.cores=4 )

lapply(pelvic.sgrids, function (x) { fcont( x, fixed="ks", col='red', add=TRUE, lty=2 ) } )
lapply(rib.sgrids, function (x) { fcont( x, fixed="ks", col='blue', add=TRUE, lty=2 ) } )
lapply(sub.pelvic.sgrids, function (x) { fcont( x, fixed="ks", col='green', add=TRUE, lty=2 ) } )


###
# MCMC

havedata <- have.pelvic
datavec <- pelvic.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[1:3]) ) )
pelvic.mcrun <- metrop( lud, initial=initpar[thesepars], nbatch=100, blen=1, scale=.1 )
pelvic.mcrun <- metrop( pelvic.mcrun, nbatch=10000, blen=1, scale=.1 )

havedata <- have.rib
datavec <- rib.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
rib.mcrun <- metrop( lud, initial=initpar[thesepars], nbatch=100, blen=1, scale=.1 )
rib.mcrun <- metrop( rib.mcrun, nbatch=10000, blen=1, scale=.1 )

# restrict to complete observations

havedata <- have.both
datavec <- pelvic.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
sub.pelvic.mcrun <- metrop( lud, initial=initpar[thesepars], nbatch=100, blen=1, scale=.1 )
sub.pelvic.mcrun <- metrop( sub.pelvic.mcrun, nbatch=10000, blen=1, scale=.1 )

havedata <- have.both
datavec <- rib.speciesdiff[ut][havedata]
stopifnot( is.finite( lud(initpar[thesepars]) ) )
sub.rib.mcrun <- metrop( lud, initial=initpar[thesepars], nbatch=100, blen=1, scale=.1 )
sub.rib.mcrun <- metrop( sub.rib.mcrun, nbatch=10000, blen=1, scale=.1 )

save( pelvic.mcrun, rib.mcrun, sub.pelvic.mcrun, sub.rib.mcrun, lud, make.spmat, pelvic.speciesdiff, rib.speciesdiff, file="shape-mcmc-results.RData" )

# all the data
layout(1:3)
for (k in 1:3) {
    tmp <- hist(c(pelvic.mcrun$batch[,k],rib.mcrun$batch[,k]),breaks=50,plot=FALSE)
    hist(rib.mcrun$batch[,k],breaks=50,col=adjustcolor('blue',.5),main=names(initpar)[k],xlim=range(tmp$breaks))
    hist(pelvic.mcrun$batch[,k],breaks=50,col=adjustcolor('red',.5),add=TRUE)
}

# complete obs
layout(1:3)
for (k in 1:3) {
    tmp <- hist(c(sub.pelvic.mcrun$batch[,k],sub.rib.mcrun$batch[,k]),breaks=50,plot=FALSE)
    hist(sub.rib.mcrun$batch[,k],breaks=50,col=adjustcolor('blue',.5),main=names(initpar)[k],xlim=range(tmp$breaks))
    hist(sub.pelvic.mcrun$batch[,k],breaks=50,col=adjustcolor('red',.5),add=TRUE)
}

if (FALSE) {

    ####
    # update branch lengths function
    # internal branches setup
    internal.lengths <- tree$edge.length
    internal.lengths[ tip.edges ] <- 0
    # and, the tips
    tip.lengths <- tree$edge.length
    tip.lengths[ - tip.edges ] <- 0

    # given parameters get rescaled branch lengths
    scale.brlens <- function (par) {
        # par = sigma2S, gammaP, xi2P
        sigma2S <- par[1]
        gammaP <- par[2]
        xi2P <- par[3]
        return( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
    }

    treemat <- treedist( adjtree, edge.length=scale.brlens(initpar[2:4]), descendants=descendants )

    # moment estimate
    have.dpelvic <- !is.na(pelvicdiff)
    have.both <- !is.na(pelvicdiff) & !is.na(ribdiff)
    # usethese <- have.both
    usethese <- have.dpelvic
    datavec <- pelvicdiff[usethese] / initpar[1]
    f1 <- function (spar) {
        sigma2S <- spar[1]
        gammaP <- spar[2]
        xi2P <- spar[3]
        elens <- ( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
        treemat <- treedist( adjtree, edge.length=elens, descendants=descendants )
        return( sum( ( datavec - treemat[usethese] )^2 ) )
    }

    est.vals <- optim( par=initpar[2:4], fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=100) )
    est.vals <- optim( par=est.vals$par, fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=10) )
    est.vals <- optim( par=est.vals$par, fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=1000) )
    # Using everything:
    # only took 31 iterations
    #    sigma2S     gammaP       xi2P 
    # 0.11451506 0.02237210 0.03384165 
    #
    # Using only ones with ribs also:
    #    sigma2S     gammaP       xi2P 
    # 0.09869252 0.01409928 0.03361052 

    new.treemat <- treedist( adjtree, edge.length=scale.brlens(est.vals$par), descendants=descendants )

    if (interactive()) {
        layout(t(1:2))
        plot( treemat[usethese], datavec )
        abline(0,1)
        plot( new.treemat[usethese], datavec )
        abline(0,1)
    }
}

if (FALSE) {

##############
# pseudolikelihood NOT WORKING
have.dpelvic <- !is.na(pelvicdiff)
chisq.datavec <- pelvicdiff[have.dpelvic]
f2 <- function (spar) {
    sigma2S <- spar[1]
    gammaP <- spar[2]
    xi2P <- spar[3]
    elens <- ( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
    treemat <- treedist( adjtree, edge.length=elens, descendants=descendants )
    return( (-1)*sum( dchisq( chisq.datavec/treemat[have.dpelvic], initpar[1], log=TRUE ) ) )
}

est.vals <- optim( par=initpar[2:4], fn=f2, method="BFGS", control=list(fnscale=1e4,maxit=10) )
est.vals <- optim( par=est.vals$par, fn=f2, method="BFGS", control=list(fnscale=1e4,maxit=10) )

}
