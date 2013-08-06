#!/usr/bin/R

#####
# MCMC
require(Matrix)

load("mcmc-setup.RData")
load("thedata-and-covmatrices.Rdata")
havedata <- !is.na(thedata)
datavec <- normdata

require(mcmc)

# return positive log-likelihood times posterior
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, delta
# priors on these are exponential
prior.means <- c(3,3,3,3,3,.2,.2,.2,.2,.2,1)
lud <- function (par) {
    if (any(par<=0)) { return( -Inf ) }
    fullmat <- make.fullmat( par )[havedata,havedata]
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( (-1) * sum( par * prior.means ) - sum( backsolve( fchol, datavec )^2 )/2 - sum(log(diag(fchol))) ) 
}

run.id <- sample.int(9999,size=1)
set.seed(run.id)

if (interactive()) {
    mcrun <- metrop( lud, initial=initpar, nbatch=10, blen=1, scale=prior.means/30 )
} else {
    mcrun <- metrop( lud, initial=initpar, nbatch=6000, blen=1, scale=prior.means/30 )
}

save(mcrun, run.id, prior.means, lud, file=paste("mcmc-run-",run.id,".RData",sep=''))
