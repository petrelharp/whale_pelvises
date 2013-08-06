#!/usr/bin/R

require(Matrix)

# treedist <- as.matrix( read.csv( file="all-sample-treedist.csv", header=TRUE) )
# tipdist <- as.matrix( read.csv( file="all-sample-tipdist.csv", header=TRUE) )
# thedata <- as.matrix( read.csv( file="all-data-rejiggered.csv", header=TRUE, row.names=1 ) )
load("thedata-and-covmatrices.Rdata")
havedata <- !is.na(thedata)
stopifnot( all( is.na(thedata) == is.na(normdata) ) )
## we only really need this components of (I-W):
## pmat <- projmatrix[1:n.tree.tips,1:n.tree.tips]
# ... but leave well enough along:
pmat <- projmatrix[havedata,]

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
    # will want to use 
    # submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
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
stopifnot( all( eigen( fullmat[havedata,havedata] )$values > -1e-8 ) )
submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
fchol <- chol( submat )

datavec <- normdata
# return negative log-likelihood for gaussian:
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, delta
llfun <- function (par) {
    fullmat <- make.fullmat( par )[havedata,havedata]
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( sum( backsolve( fchol, datavec )^2 )/2 + sum(log(diag(fchol))) ) 
}
stopifnot( is.finite(llfun(initpar)) )

save(make.fullmat, initpar, thedata, normdata, species.Pmat, species.Qmat, sample.Pmat, sample.Qmat, sample.Pcoef, species.transmat, sample.transmat, file="mcmc-setup.RData")

if (!file.exists("analysis-results.RData")) {
    do.parallel <- TRUE
    if (do.parallel) {
        require(parallel)
        pjobs <-  list( mcparallel( { optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( maxit=1000 ) ) } ),
                mcparallel( { optim( par=initpar, fn=llfun, method="BFGS", control=list( maxit=200 ) ) } ),
                mcparallel( { optim( par=initpar, fn=llfun, method="L-BFGS-B", control=list( trace=3 ), lower=1e-3 ) } )
            )
        mlestims <- mccollect( pjobs, wait=TRUE )
        mlestim1 <- mlestims[[1]]
        mlestim2 <- mlestims[[2]]
        mlestim3 <- mlestims[[3]]
    } else {
        mlestim1 <- optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( maxit=1000 ) )
        mlestim2 <- optim( par=initpar, fn=llfun, method="BFGS", control=list( maxit=200 ) )
        mlestim3 <- optim( par=initpar, fn=llfun, method="L-BFGS-B", control=list( ), lower=1e-3 )
        mlestims <- list( mlestim1, mlestim2, mlestim3 )
    }

    save( mlestims, file="analysis-results.RData" )
}

if (FALSE) {
    ####
    # examine results
    load("analysis-results.RData")

    mlpars <- as.data.frame( rbind( initpar, do.call( rbind, lapply(mlestims,"[[","par") ) ) )
    mlpars$ll <- apply( mlpars, 1, llfun )

    mlfullmat <- make.fullmat(mlestim1$par)

    source("correlated-traits-fns.R")

    ##
    # look at residuals...
    no.lefties <- thedata
    no.lefties[,grep("left.",colnames(thedata))] <- NA
    pred.leftbones <- predgaus( as.vector(no.lefties), unlist(phylomeans)[col(thedata)], mlfullmat )
    dim(pred.leftbones) <- dim(thedata)
    resid.leftbones <- (thedata-pred.leftbones)[,grep("left.",colnames(thedata))]

    ### 
    # diff betweens individual's ribs:
    boneinds <- grep(".pelvic",colnames(thedata),fixed=TRUE)
    bonediffmat <- sapply( levels(whales$specimen), function (sp) {
                side <- matrix(0,nrow=nrow(thedata),ncol=ncol(thedata))
                side[col(thedata) == boneinds[1]] <- 1
                side[col(thedata) == boneinds[2]] <- -1
                out <- as.numeric( ( rownames(thedata)[row(thedata)] == sp ) ) * side
                if( any( is.na( thedata[ out!=0 ] ) ) ) { 
                    return(NULL) 
                } else {
                    return( out )
                }
            } )
    bonediffmat <- do.call(rbind,lapply(bonediffmat,as.vector))
    zdata <- thedata
    zdata[is.na(zdata)] <- 0
    bonediffs <- bonediffmat %*% as.vector( zdata )
    summary(bonediffs)
    sqrt(var(bonediffs))


    ml.bonediff.covmat <- tcrossprod( bonediffmat %*% mlfullmat, bonediffmat )
    range( ml.bonediff.covmat[ row(ml.bonediff.covmat) > col(ml.bonediff.covmat) ] )
    range( diag(ml.bonediff.covmat) )

    bonediff.covmat <- tcrossprod( bonediffmat %*% fullmat, bonediffmat )
    range( bonediff.covmat[ row(bonediff.covmat) > col(bonediff.covmat) ] )
    range( diag(bonediff.covmat) )



    subsp <- c("STENELLA_LONGIROSTRIS","STENELLA_FRONTALIS","STENELLA_COERULEOALBA","DELPHINUS_DELPHIS")
    subind <- c( match( subsp, rownames(thedata) ), match( whales$specimen[ whales$species %in% subsp ], rownames(thedata) ) )
    submatind <- as.vector( outer( subind, (0:(ncol(thedata)-1))*nrow(thedata), "+" ) )

    image( Matrix( cov2cor( mlfullmat ) ) )
    image( Matrix( cov2cor( mlfullmat[submatind,submatind] ) ) )

    image( Matrix( cov2cor( fullmat ) ) )
    image( Matrix( cov2cor( fullmat[submatind,submatind] ) ) )

    layout(matrix(1:36,nrow=6))
    par(mar=c(0,0,0,0)+.1)
    for (i in 0:5) for (j in 0:5) {
        i.thisind <- submatind[ (submatind-1)%%6 == i ]
        j.thisind <- submatind[ (submatind-1)%%6 == j ]
        image( mlfullmat[i.thisind,j.thisind] )
    }


    mlchol <- chol(mlfullmat[havedata,havedata])
    normresids <- backsolve( mlchol, datavec )
    layout(t(1:2))
    plot( normresids )
    outliers <- identify( normresids )
    outliercoef <- mlchol[outliers,]

    plot( thedata[havedata] )
    identify( thedata[havedata], labels=paste(rownames(thedata)[row(thedata)],colnames(thedata)[col(thedata)],sep='.')[havedata] )
}

#### TO-DO:
## SUBTRACT OFF MEAN TRAIT VALUES
## NORMALIZE BY SEXUAL DIMORPHISM
## REVISIT SAME DELTA FOR TESTES, RIB, and PELVIS

####
# testing:


if (FALSE) {  # inference simulating under the model works:
    require(mvtnorm)
    fakedata <- fchol %*% rnorm(nrow(fchol))
    fake.llfun <- function (par) {
        fchol <- chol(make.fullmat(par)[havedata,havedata])
        return( sum( backsolve( fchol, fakedata )^2 )/2 + sum(log(diag(fchol)))/2 ) 
    }
    fake.mlestim <- optim( par=initpar, fn=fake.llfun, method="L-BFGS-B", control=list( fnscale=1e3, trace=3 ), lower=1e-3 )
    rbind( initpar, fake.mlestim$par ) # looks good
}
