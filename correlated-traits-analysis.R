#!/usr/bin/R

require(Matrix)


#####
# set up the model

# associate matrices with each internal branch of the tree:
#  these parameters are: theta = sigmaL, betaT, betaP, sigmaR, sigmaP.
# the variables are, in order: Length, Testes, Rib-left, Rib-right, Pelvis-left, Pelvis-right
# Also, associate matrices with tip-edges.
#  the paramters are: zetaL, zetaR, omegaR, zetaP, omegaP
species.paramnames <- matrix(
              c('sigmaL','','','','','',
                'sigmaLdeltaT','betaT','','','','',
                'sigmaLdeltaP','betaR','sigmaR','','','',
                'sigmaLdeltaP','betaR','sigmaR','','','',
                'sigmaLdeltaR','betaP','','sigmaP','','',
                'sigmaLdeltaR','betaP','','sigmaP','',''), nrow=6, byrow=TRUE )
sample.paramnames <- matrix( 
              c('zetaL','','','','','',
                '','','','','','',
                'zetaLdeltaP','','zetaR','-omegaR','','',
                'zetaLdeltaP','','zetaR','omegaR','','',
                'zetaLdeltaR','','','','zetaP','-omegaP',
                'zetaLdeltaR','','','','zetaP','omegaP'), nrow=6, byrow=TRUE )
species.justparams <- gsub("\\..*","",gsub("^[-+]","",species.paramnames))
species.deltas <- gsub("[^.]\\.","",species.paramnames)
species.varnames <- setdiff( unique( species.justparams ), '' )
species.Tind <- which( species.justparams != '' )  # nonzero elements of species.transmat
species.Sind <- match( species.justparams[species.Tind], species.varnames )     # put params[species.Sind] into species.transmat[species.Tind]
stopifnot( length(species.Sind)==length(species.Tind) )
sample.justparams <- gsub("\\..*","",gsub("^[-+]","",sample.paramnames))
sample.signs <- ifelse( grepl( "^-", sample.paramnames ), -1, +1 )
sample.varnames <- setdiff( unique( sample.justparams ), '' )
sample.Tind <- which( sample.justparams != '' )
sample.Sind <- match( sample.justparams[sample.Tind], sample.varnames )     # put params[sample.Sind]*sample.Pcoef into sample.transmat[sample.Tind]
sample.Pcoef <- sample.signs[ sample.Tind ]
stopifnot( length(sample.Sind)==length(sample.Tind) )


### initial values
# from initial-values.R
initpar <- c(
        sigmaL=3.16,
        betaT=6.5,
        betaP=2,
        betaR=.1,
        sigmaR=.5,
        sigmaP=1.1,
        zetaL=.05,
        zetaR=.06,
        omegaR=.01,
        zetaP=.12,
        omegaP=.02,
        # deltaT=sqrt(1.4),
        # deltaP=sqrt(1.3),
        # deltaR=sqrt(1.2),
        sigmaLdeltaT=3.16*sqrt(1.4),
        sigmaLdeltaP=3.16*sqrt(1.3),
        sigmaLdeltaR=3.16*sqrt(1.2),
        zetaLdeltaP=.05*sqrt(1.3),
        zetaLdeltaR=.05*sqrt(1.2)
    )
stopifnot( all(names(initpar) %in% c( species.justparams, sample.justparams ) ) )

species.parindices <- match( species.varnames, names(initpar) )
sample.parindices <- match( sample.varnames, names(initpar) )
species.params <- initpar[species.parindices]
sample.params <- initpar[sample.parindices]
# delta <- initpar[11:13]
species.transmat <- matrix(0,nrow=nrow(species.paramnames),ncol=ncol(species.paramnames))
sample.transmat <- matrix(0,nrow=nrow(sample.paramnames),ncol=ncol(sample.paramnames))
species.transmat[species.Tind] <- species.params[species.Sind] # * c(1,delta)[1+species.Dind]
sample.transmat[sample.Tind] <- sample.params[sample.Sind] * sample.Pcoef #*c(1,delta)[1+sample.Dind]
species.covmat <- as.matrix( tcrossprod(species.transmat) )
sample.covmat <- as.matrix( tcrossprod(sample.transmat) )
sptransmat <- species.transmat
samtransmat <- sample.transmat


#####
# Read in the data

load("thedata-and-covmatrices.Rdata")
havedata <- !is.na(thedata)
pmat <- projmatrix[havedata,]

###
# function to construct full covariance matrix
make.fullmat <- function (par) {
    # return full covariance matrix for all data (observed and unobserved)
    species.params <- par[species.parindices]
    sample.params <- par[sample.parindices]
    # delta <- par[10+1:3]
    species.transmat[species.Tind] <- species.params[species.Sind] #* c(1,delta)[1+species.Dind]
    sample.transmat[sample.Tind] <- sample.params[sample.Sind] * sample.Pcoef #* c(1,delta)[1+sample.Dind]
    species.covmat <- as.matrix( tcrossprod(species.transmat) )
    sample.covmat <- as.matrix( tcrossprod(sample.transmat) )
    fullmat <-  kronecker( species.covmat, species.treemat ) + kronecker( sample.covmat, sample.treemat )
    # will want to use 
    # submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    return( fullmat )
}

# construct full matrix
fullmat <- make.fullmat( initpar )
colnames( fullmat ) <- rownames( fullmat ) <- outer( rownames(thedata), colnames(thedata), paste, sep='.' )
stopifnot( all( eigen( fullmat[havedata,havedata] )$values > -1e-8 ) )
submat <- ( ( crossprod( pmat, fullmat[havedata,havedata]) %*% pmat ) )
fchol <- chol( submat )

datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data
# return negative log-likelihood for gaussian:
#  parameters are: sigmaL, betaT, betaP, sigmaR, sigmaP, zetaL, zetaR, omegaR, zetaP, omegaP, delta
llfun <- function (par) {
    fullmat <- make.fullmat( par )[havedata,havedata]
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( sum( backsolve( fchol, datavec, transpose=TRUE )^2 )/2 + sum(log(diag(fchol))) ) 
}
stopifnot( is.finite(llfun(initpar)) )

save(make.fullmat, initpar, thedata, pmat, havedata, llfun, species.parindices, sample.parindices, species.transmat, species.Tind, species.Sind, sample.transmat, sample.Tind, sample.Sind, sample.Pcoef, species.treemat, sample.treemat, file="mcmc-setup.RData")

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
    mlpars$convergence <- c( NA, sapply( mlestims, "[[", "convergence" ) )

    mlfullmat <- make.fullmat(mlestim1$par)

    source("correlated-traits-fns.R")

    ##
    # differences versus covariance
    diffmat <- outer( as.vector(thedata), as.vector(thedata), "-" )
    layout( matrix(1:6,nrow=2) )
    for (k in 1:6) { 
        theseones <- 1:nrow(thedata) + (k-1)*nrow(thedata)
        plot( cov2cor(fullmat)[ theseones, theseones ], diffmat[ theseones, theseones ] )
    }
    load("all-sample-tree.RData")
    edge.indices <- tree$edge[,2] # associate each edge with the downstream node
    tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
    treedist <- treedist(tree,edge.length=ifelse(1:Nedge(tree) %in% tip.edges, 0, tree$edge.length))
    layout( matrix(1:6,nrow=2) )
    for (k in 1:6) { 
        theseones <- 1:nrow(thedata) + (k-1)*nrow(thedata)
        plot( treedist, cov2cor(fullmat)[ theseones, theseones ] )
    }


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
    normresids <- backsolve( mlchol, datavec, transpose=TRUE )
    layout(t(1:2))
    plot( normresids )
    outliers <- identify( normresids )
    outliercoef <- mlchol[outliers,]

    plot( thedata[havedata] )
    identify( thedata[havedata], labels=paste(rownames(thedata)[row(thedata)],colnames(thedata)[col(thedata)],sep='.')[havedata] )
}
