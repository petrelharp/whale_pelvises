#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file from previous MCMC run." ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
if (interactive()) {
    nbatches <- 10
}

new.mcmc <- is.null(infile)
old.run.id <- NA
if (!new.mcmc) { 
    load(infile) 
    old.run.id <- run.id
}

#####
# MCMC
require(Matrix)

load("mcmc-setup.RData")
load("thedata-and-covmatrices.Rdata")
havedata <- !is.na(thedata)
datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data

require(mcmc)

# return positive log-likelihood times posterior
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, sigmaLdeltaT, sigmaLdeltaP, sigmaLdeltaR, zetaLdeltaP, zetaLdeltaR, 
# priors on these are exponential
prior.means <- c(3,3,3,3,3,.1,.1,.1,.1,.1,1,1,1,.1,.1)
stopifnot( length(prior.means) == length(initpar) )
lud <- function (par) {
    if (any(par<=0)) { return( -Inf ) }
    fullmat <- make.fullmat( par )[havedata,havedata]
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( (-1) * sum( par * prior.means ) - sum( backsolve( fchol, datavec, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}

run.id <- sample.int(9999,size=1)
set.seed(run.id)

if (new.mcmc) {
    mcrun <- metrop( lud, initial=initpar, nbatch=nbatches, blen=1, scale=prior.means/100 )
} else {
    mcrun <- metrop( mcrun, nbatch=nbatches, blen=1, scale=prior.means/100 )
}

save(mcrun, run.id, old.run.id, prior.means, lud, file=paste("mcmcs/mcmc-run-",run.id,".RData",sep=''))
