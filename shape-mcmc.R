#/usr/bin/R --vanilla
require(optparse)

usage <- "MCMC over parameters in shape analysis"
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
basedir <- paste(type,"shape-mcmc",sep='-')
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
initpar <- c( kS=20, sigma2S=2.5, gammaP=0.5, xi2P=0.03 )
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
    return( (-1) * sum( (par / prior.means) ) - sum( backsolve( fchol, x, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}
# priors: exponential
prior.means <- c(20,2,2)
stopifnot( length(prior.means) == length(thesepars) )

have.pelvic <- !is.na(pelvic.speciesdiff[ut])
have.rib <- !is.na(rib.speciesdiff[ut])
have.both <- have.pelvic & have.rib
havedata <- have.both

###
# MCMC

if (type=='pelvic') {
    havedata <- have.pelvic
    datavec <- pelvic.speciesdiff[ut][havedata]
} else if (type == "rib") {
    havedata <- have.rib
    datavec <- rib.speciesdiff[ut][havedata]
} else if (type == "sub.pelvic") {
    havedata <- have.both
    datavec <- pelvic.speciesdiff[ut][havedata]
}

stopifnot( is.finite( lud(initpar[1:3]) ) )

if (new.mcmc) {
    mcrun <- metrop( lud, initial=initpar[thesepars], nbatch=nbatches, blen=blen, scale=.1 )
} else {
    mcrun <- metrop( mcrun, nbatch=nbatches, blen=blen, scale=.1 )
}
colnames(mcrun$batch) <- thesepars


save( mcrun, lud, havedata, datavec, make.spmat, opt, file=outfile )
