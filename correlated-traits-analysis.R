#!/usr/bin/R

require(Matrix)

# treedist <- as.matrix( read.csv( file="all-sample-treedist.csv", header=TRUE) )
# tipdist <- as.matrix( read.csv( file="all-sample-tipdist.csv", header=TRUE) )
# thedata <- as.matrix( read.csv( file="all-data-rejiggered.csv", header=TRUE, row.names=1 ) )
load("thedata-and-covmatrices.Rdata")
havedata <- !is.na(thedata)

# associate P and Q with each internal branch of the tree:
#  these parameters are: theta = sigmaL, betaT, betaP, sigmaR, sigmaP.
#  The vector Pmat is such that theta[Pmat] gives the NONZERO elements of the corresponding sqrt-covariance matrix.
#  ... which are (branch length) * (sigmaL, sigmaL, sigmaL, sigmaL, betaT, betaP, sigmaR, sigmaP)
#  and the matrix Q is similar, except says where to put delta (on all sigmaL but the first one)
species.Pmat <- c(1,1,1,1,1,1,2,3,3,4,4,5,5)
species.Qmat <- c(0,1,1,1,1,1,0,0,0,0,0,0,0)

# Now add tips for samples:
# When constructing P, Q, 
#  the paramters are: zetaL, zetaR, omegaR, zetaP, omegaP
# and the nonzero elements, in order, are
#  (zetaL, zetaL, zetaL, zetaL, zetaL, zetaR, zetaR, omegaR, -omegaR, zetaP, zetaP, omegaP, -omegaP)
# again, put delta on all zetaL but the first one
#  now, nonzero elements are theta[Pmat] * Pcoef
sample.Pmat <- c(1,1,1,1,1,2,2,3,3,4,4,5,5)
sample.Pcoef <- c(1,1,1,1,1,1,1,1,-1,1,1,1,-1)
sample.Qmat <- c(0,1,1,1,1,0,0,0,0,0,0,0,0)

# set up with some initial parameters
# the variables are, in order: Length, Testes, Rib-left, Rib-right, Pelvis-left, Pelvis-right
species.transmat <- Matrix( 
              c(1,0,0,0,
                1,1,0,0,
                1,0,1,0,
                1,0,1,0,
                1,1,0,1,
                1,1,0,1), nrow=6, byrow=TRUE, sparse=TRUE )
sample.transmat <- Matrix( 
              c(1,0,0,0,0,
                0,0,0,0,0,
                1,1,1,0,0,
                1,1,-1,0,0,
                1,0,0,1,1,
                1,0,0,1,-1), nrow=6, byrow=TRUE, sparse=TRUE )

make.fullmat <- function (par) {
    # return full covariance matrix for all data (observed and unobserved)
    #  parameters are: ( sigmaL, betaT, betaP, sigmaR, sigmaP ), (zetaL, zetaR, omegaR, zetaP, omegaP), (delta)
    species.params <- par[1:5]
    sample.params <- par[5+1:5]
    delta <- par[11]
    species.transmat@x <- as.vector( ( species.params[species.Pmat] ) * ( 1 + species.Qmat * ((delta)-1) ) )
    sample.transmat@x <- as.vector( ( sample.params[sample.Pmat] * sample.Pcoef ) * ( 1 + sample.Qmat * ((delta)-1) ) )
    species.covmat <- as.matrix( tcrossprod(species.transmat) )
    sample.covmat <- as.matrix( tcrossprod(sample.transmat) )
    fullmat <-  kronecker( species.covmat, species.treemat ) + kronecker( sample.covmat, sample.treemat )
    return( fullmat )
}

### initial values
# from initial-values.R
initpar <- c(
        sigmaL=3.16,
        betaT=6.5,
        betaP=2,
        sigmaR=.5,
        sigmaP=1.1,
        zetaL=.05,
        zetaR=.06,
        omegaR=.01,
        zetaP=.12,
        omegaP=.02,
        delta=sqrt(1.3)
    )
# construct full matrix
fullmat <- make.fullmat( initpar )
colnames( fullmat ) <- rownames( fullmat ) <- outer( rownames(thedata), colnames(thedata), paste, sep='.' )
fchol <- chol( fullmat[ havedata, havedata ] )

datavec <- thedata[havedata]
norm.datavec <- sweep( thedata, 2, colMeans(thedata,na.rm=TRUE), "-" )[havedata]
# return negative log-likelihood for gaussian:
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, delta
llfun <- function (par) {
    fchol <- chol(make.fullmat(par)[havedata,havedata])
    return( sum( backsolve( fchol, norm.datavec )^2 )/2 + sum(log(diag(fchol)))/2 ) 
}

mlestim1 <- optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( fnscale=1e7, trace=3 ) )
mlestim2<- optim( par=initpar, fn=llfun, method="BFGS", control=list( fnscale=1e7, trace=3 ) )
mlestim3 <- optim( par=initpar, fn=llfun, method="L-BFGS-B", control=list( fnscale=1e7, trace=3 ), lower=1e-3 )

mlestims <- data.frame( rbind(initpar, mlestim1$par, mlestim2$par, mlestim3$par ) )
mlestims$ll <- apply( mlestims, 1, llfun )

require(mcmc)
# return positive log-likelihood times posterior
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, delta
# priors on these are exponential
prior.means <- c(3,3,3,3,3,.2,.2,.2,.2,.2,1)
lud <- function (par) {
    if (any(par<=0)) { return( -Inf ) }
    fchol <- chol(make.fullmat(par)[havedata,havedata])
    return( (-1) * sum( par * prior.means ) - sum( backsolve( fchol, norm.datavec )^2 )/2 - sum(log(diag(fchol)))/2 ) 
}

mcrun <- metrop( lud, initial=initpar, nbatch=100, blen=1, scale=prior.means/10 )


#### TO-DO:
## SUBTRACT OFF MEAN TRAIT VALUES
## NORMALIZE BY SEXUAL DIMORPHISM
## REVISIT SAME DELTA FOR TESTES, RIB, and PELVIS

####
# testing:

untf.data <- backsolve( fchol, norm.datavec )
badones <- which(abs(untf.data)>800)
fchol.badones <- backsolve( fchol, sapply( badones, function(k) as.numeric( ((1:nrow(fchol))==k) ) ) )
lapply(1:ncol(fchol.badones), function(k) {
        these <- which( abs(fchol.badones[,k]) > .1 )
        x <- which( havedata, arr.ind=TRUE )[these,]
        dim(x) <- c(length(x)/2,2)
        cbind( rownames(thedata)[x[,1]], colnames(thedata)[x[,2]] )
    } )

if (FALSE) {  # inference simulating under the model works:
    require(mvtnorm)
    fakedata <- fchol %*% rnorm(nrow(fchol))
    fake.llfun <- function (par) {
        fchol <- chol(make.fullmat(par)[havedata,havedata])
        return( sum( backsolve( fchol, fakedata )^2 )/2 + sum(log(diag(fchol)))/2 ) 
    }
    fake.mlestim <- optim( par=initpar, fn=fake.llfun, method="L-BFGS-B", control=list( fnscale=1e7, trace=3 ), lower=1e-3 )
    rbind( initpar, fake.mlestim$par ) # looks good
}
