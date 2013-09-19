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
if (interactive()) { type <- 'pelvic'; nbatches <- 100; blen <- 1; restart <- FALSE }
if (is.null(opt$type)) { stop('type should be one of (pelvic, rib, or sub.pelvic)') }
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

load("shape-brlens-stuff.RData")
load("spmapping.RData")
load("shared-num-bones.RData")
new.mcmc <- is.null(opt$mcmcfile)
if (!new.mcmc) { load(mcmcfile) }

#####
# Likelihood, species tree

# Estimation:
stype <- gsub("sub.","",type)
initpar <- c( initpar[[stype]][c('ks','xi2P','xi2P')]/c(1,2,2), ( 5 * sptree$edge.length ) )
parscale <- c(.5,initpar[-1]/50)
stopifnot( class(shared.paths) == "dsyMatrix"  & shared.paths@uplo == "U" )  # if so, changing upper tri also changes lower tri
ut <- upper.tri(pelvic.speciesdiffsq)
make.spmat <- function ( edgelens ) {
    shared.paths@x[ sput ][spmap.nonz] <- as.vector( sp.mapping %*% (edgelens) )  # note: update @x rather than entries to preserve symmetry
    return( shared.paths )
}

# priors: exponential
prior.means <- c(20,rep(.1,length=2+length(sptree$edge.length)))

lud <- function (par, ...) {
    # par = ks, xi2i, xi2s, (edgelens)
    ks <- par[1]
    xi2i <- par[2]
    xi2s <- par[3]
    edgelens <- par[-(1:3)]
    if (any(par<=0)) { return( -Inf ) }
    #   covariances between the means of all species-species comparisons is,
    #     with y = shared.indivs * xi2i
    #          z = shared.samples * xi2s
    #     and y.z = shared.samples.indivs * xi2i * xi2s
    #     and y.sq = shared.indivs.sq * xi2i^2
    #          z.sq = shared.samples.sq * xi2s^2
    #     shared.paths^2 + 2 * shared.paths * (y+z) + 2 * y.z + y.sq + z.sq
    spmat <- make.spmat( edgelens, ... )[usedata,usedata]
    covmat <- ( spmat^2 + 
            2 * spmat * ( shared.indivs * xi2i + shared.samples * xi2s ) + 
            2 * shared.samples.indivs * xi2i * xi2s +
            shared.indivs.sq * xi2i^2 +
            shared.samples.sq * xi2s^2 )
    x <- ( datavec - ks * diag(spmat) )
    fchol <- chol(2*ks*covmat)
    return( sum( (par / prior.means) ) - sum( backsolve( fchol, x, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}

have.pelvic <- !is.na(pelvic.speciesdiffsq[ut])
have.rib <- !is.na(rib.speciesdiffsq[ut])
have.both <- have.pelvic & have.rib

sharenames <- c( "shared.indivs", "shared.samples", "shared.samples.indivs", "shared.samples.sq", "shared.indivs.sq" )
if (type=='pelvic') {
    usedata <- have.pelvic
    datavec <- pelvic.speciesdiffsq[ut][usedata]
    for (x in sharenames) { assign( x, get(paste(type,x,sep='.'))[usedata,usedata] ) }
} else if (type == "rib") {
    usedata <- have.rib
    datavec <- rib.speciesdiffsq[ut][usedata]
    for (x in sharenames) { assign( x, get(paste(type,x,sep='.'))[usedata,usedata] ) }
} else if (type == "sub.pelvic") {
    usedata <- have.both
    datavec <- pelvic.speciesdiffsq[ut][usedata]
    for (x in sharenames) { assign( x, get(paste(type,x,sep='.'))[usedata,usedata] ) }
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


save( mcrun, lud, usedata, datavec, make.spmat, opt, file=outfile )
