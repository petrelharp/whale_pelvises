#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default='', help=".RData file from previous MCMC run." ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-o","--outdir"), type="character", default='mcmcs/', help="Directory to put output in. [default \"%default\"]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
if (interactive()) {
    nbatches <- 10
}

new.mcmc <- (infile == '')
if (!new.mcmc) { load(infile) }
old.run.id <- if (new.mcmc) { NA } else { run.id }

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
# priors are Gaussian with these SDs
prior.means <- c(3,3,3,3,3,3,.1,.1,.1,.1,.1,1,1,1,.1,.1)
# constrain these to be nonnegative:
nonnegs <- c("sigmaL", "betaT", "sigmaR", "sigmaP", "zetaL", "zetaR", "omegaR", "zetaP", "omegaP" )
nonneg.inds <- match( nonnegs, names(initpar) ) 
stopifnot( length(prior.means) == length(initpar) )
lud <- function (par) {
    if (any(par[nonneg.inds]<=0)) { return( -Inf ) }
    fullmat <- make.fullmat( par )[havedata,havedata]
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( (-1/2) * sum( (par / prior.means)^2 ) - sum( backsolve( fchol, datavec, transpose=TRUE )^2 )/2 - sum(log(diag(fchol))) ) 
}

run.id <- sample.int(9999,size=1)
set.seed(run.id)

par.scale <- prior.means/100
par.scale[ names(initpar) %in% c( "sigmaLdeltaT", "sigmaLdeltaP", "sigmaLdeltaR" ) ] <- 1/10

if (new.mcmc) {
    mcrun <- metrop( lud, initial=initpar, nbatch=nbatches, blen=1, scale=prior.means/100 )
} else {
    mcrun <- metrop( mcrun, nbatch=nbatches, blen=1, scale=prior.means/100 )
}

save(mcrun, run.id, old.run.id, prior.means, lud, file=paste(outdir, "mcmc-run-",run.id,".RData",sep=''))
