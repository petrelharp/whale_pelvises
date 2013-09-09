#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(Matrix)

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
    spmat <- make.spmat( par[2:3], ... )[usedata,usedata]
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
usedata <- have.both

###
# fix kS and look at 2D density
#  across parameter grid and sampled values in sampled.edge.testes

ks.vals <- seq(17,20,length.out=4)
sigma2S.vals <- seq(1,5,length.out=40)
gammaP.vals <- seq(-1,2,length.out=40)
pargrid <- expand.grid( ks=ks.vals, sigma2S=sigma2S.vals, gammaP=gammaP.vals )

if (interactive()) {
    usedata <- have.pelvic
    datavec <- pelvic.speciesdiff[ut][usedata]
    stopifnot( is.finite( lud(initpar[thesepars]) ) )
    pargrid$pelvic <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

    usedata <- have.ribs
    datavec <- rib.speciesdiff[ut][usedata]
    stopifnot( is.finite( lud(initpar[thesepars]) ) )
    pargrid$ribs <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

    usedata <- ( have.pelvic & have.ribs )
    datavec <- pelvic.speciesdiff[ut][usedata]
    stopifnot( is.finite( lud(initpar[thesepars]) ) )
    pargrid$sub.pelvic <- apply( pargrid, 1, function (x) lud( x[1:3] ) )

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

    fcont("pelvic", fixed="sigma2S",xlab='kS',ylab='gammaP')
    fcont("ribs", fixed="sigma2S",add=TRUE,col='blue',xlab='kS',ylab='gammaP')
    fcont("sub.pelvic", fixed="sigma2S",add=TRUE,col='green',xlab='kS',ylab='gammaP')
}

# effect of sampling from posterior on testes size?
nreps <- 100
do.parallel <- FALSE
if (do.parallel) {
    require(parallel)
    this.lapply <- function (...) { mclapply( ..., mc.cores=8 ) }
} else {
    this.lapply <- lapply
}

sampled.edge.testes <- t( replicate(nreps, sample.values() ) )
sampled.sp.edge.testes <- sampled.edge.testes[ , tree.translate ]

usedata <- have.pelvic
datavec <- pelvic.speciesdiff[ut][usedata]
pelvic.sgrids <- this.lapply( 1:nreps, function (k) {
        # k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    } )

usedata <- have.ribs
datavec <- rib.speciesdiff[ut][usedata]
rib.sgrids <- this.lapply( 1:nreps, function (k) {
        # k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    } )

usedata <- ( have.pelvic & have.ribs )
datavec <- pelvic.speciesdiff[ut][usedata]
sub.pelvic.sgrids <- this.lapply( 1:nreps, function (k) {
        # k <- sample(1:nrow(sampled.sp.edge.testes),1)
        edge.weights <- sampled.sp.edge.testes[k,]
        sub.pelvic.sgrid <- apply( pargrid, 1, function (x) lud( x[1:3], edge.weights=edge.weights ) )
    } )

run.id <- sprintf( sample(1e4,1), fmt="%04.0f" )
if (no.females) { run.id <- paste("no-females-",run.id,sep='') }
save( pargrid, pelvic.sgrids, rib.sgrids, sub.pelvic.sgrids, make.spmat, thesepars, lud, prior.means, sampled.edge.testes, sampled.sp.edge.testes, file=paste("likelihood-surface-", run.id, ".RData",sep='') )
