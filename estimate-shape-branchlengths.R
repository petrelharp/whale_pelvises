#/usr/bin/R --vanilla
require(optparse)

usage <- "MCMC over shape branch lengths."
option_list <- list(
        make_option( c("-t","--type"), type="character", default=NULL, help="Which dataset? (pelvic, rib, or sub.pelvic)"),
        make_option( c("-m","--mcmcfile"), type="character", default=NULL, help=".RData file containing previous MCMC run." ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=1, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
if (interactive()) { type <- 'pelvic'; nbatches <- 100; blen <- 10; restart <- FALSE }
basedir <- paste(type,"shape-brlens",sep='-')
basename <- paste(basedir, sprintf( sample(1e4,1), fmt="%04.0f" ), sep='/')

outfile <- paste(basename, "-mcmc-run.RData",sep='')
if (logfile=="") { logfile <- paste(basename,"-mcmc-run.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="output", split=interactive()) 
    if (!interactive()) { sink(file=logcon, type="message",append=TRUE) }
}

source("correlated-traits-fns.R")
require(ape)
require(Matrix)
require(mcmc)

load("shape-stuff.RData")
new.mcmc <- is.null(opt$mcmcfile)
if (!new.mcmc) { load(mcmcfile) }

#####
# Likelihood, species tree

# Estimation:
initpar <- c( kS=20, xi2P=, ( 5 * sptree$edge.length ) )
parscale <- c(.5,initpar[-1]/50)
stopifnot( class(shared.paths) == "dsyMatrix"  & shared.paths@uplo == "U" )  # if so, changing upper tri also changes lower tri
make.spmat <- function ( edgelens ) {
    shared.paths@x[ sput ][spmap.nonz] <- as.vector( sp.mapping %*% (edgelens) )  # note: update @x rather than entries to preserve symmetry
    return( shared.paths )
}
lud <- function (par, ...) {
    # par = ks, (edgelens)
    ks <- par[1]
    edgelens <- par[-1]
    if (any(par<=0)) { return( -Inf ) }
    spmat <- make.spmat( edgelens, ... )[usedata,usedata]
    x <- ( datavec - ks * diag(spmat) )
    fchol <- chol(2*ks*spmat^2)
    return( sum( (par / prior.means) ) - sum( backsolve( fchol, x, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}
# priors: exponential
prior.means <- c(20,rep(.1,length=length(sptree$edge.length)))

have.pelvic <- !is.na(pelvic.speciesdiff[ut])
have.rib <- !is.na(rib.speciesdiff[ut])
have.both <- have.pelvic & have.rib


if (type=='pelvic') {
    usedata <- have.pelvic
    datavec <- pelvic.speciesdiff[ut][usedata]
} else if (type == "rib") {
    usedata <- have.rib
    datavec <- rib.speciesdiff[ut][usedata]
} else if (type == "sub.pelvic") {
    usedata <- have.both
    datavec <- pelvic.speciesdiff[ut][usedata]
} else {
    stop( "type must be 'pelvic', 'rib', or 'sub.pelvic'." )
}

stopifnot( is.finite( lud(initpar) ) )

if (new.mcmc) {
    mcrun <- metrop( lud, initial=initpar, nbatch=nbatches, blen=blen, scale=parscale )
} else {
    mcrun <- metrop( mcrun, nbatch=nbatches, blen=blen, scale=parscale )
}
colnames(mcrun$batch) <- names(initpar)


save( mcrun, lud, havedata, datavec, make.spmat, opt, file=outfile )
