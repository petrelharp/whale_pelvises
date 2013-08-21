#!/usr/bin/R

require(ape)
require(Matrix)

basedir <- scriptdir <- ".."
source(paste(scriptdir,"correlated-traits-fns.R",sep="/"))

tree_file <- paste(basedir, "consensusTree_ALL_CETACEA.tree", sep='/')
species_tree<-read.nexus(file=tree_file)
load( paste(basedir,"all-sample-tree.RData",sep='/') ) # gets tree

whales <- read.csv(paste(basedir,"whales.csv",sep='/'),header=TRUE)
allspecies <- intersect(tree$tip.label, species_tree$tip.label)

do.simple <- TRUE
if (do.simple) {
    ## SIMPLE EXAMPLE
    allspecies <- c( "PHOCOENOIDES_DALLI", "STENELLA_ATTENUATA", "STENELLA_LONGIROSTRIS" )
    whales <- droplevels( subset(whales,species%in%allspecies) )
    tree <- drop.tip( tree, setdiff( tree$tip.label, c( allspecies, as.character(whales$specimen) ) ) )

    bones <- read.table(paste(basedir,"62_add_centroids.out",sep='/'), header=TRUE)
    bones <- droplevels( subset(bones, species %in% allspecies & ! (specimen == "LACM_54109") ) )
    bones <- bones[ setdiff( colnames(bones), "absolute_volume" ) ]
    species <- read.table(paste(basedir,"52_sexual_dimorphism.out",sep='/'), header=TRUE)
    morphology <- read.table(paste(basedir,"morphology_table_2013_June_27-plr.txt",sep='/'), sep='\t', header=TRUE)
    stopifnot( all.equal( allspecies , sort( unique( levels(bones$species) ) ) ) )
    bones$species <- factor( bones$species, levels=allspecies )
    tmp <- data.frame( 
            species=allspecies,
            bodylength=tapply( bones$bodylength, bones$species, mean )
            )
    stopifnot( any( tapply( bones$bodylength, bones$species, var )>0, na.rm=TRUE ) )
    species <- subset( species, species %in% allspecies )
    species$species <- factor( species$species, levels=allspecies )
    species <- merge( species, tmp, by="species", all.x=TRUE, all.y=TRUE )
    morphology$species <- factor( toupper(morphology$species), levels=allspecies )
    species <- merge( species, morphology, by='species', all.x=TRUE )

    # get thedata ready
    thedata <- matrix( c(NA), nrow=Nnode(tree)+Ntip(tree), ncol=6 )
    colnames(thedata) <- c("bodylength","actual_testes_mass_max","left.rib","right.rib","left.pelvic","right.pelvic")
    rownames(thedata) <- c( tree$tip.label, paste("NA",1:Nnode(tree),sep='.') )
    species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
    sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
    thedata[species.nodes,intersect(names(species),colnames(thedata))] <- as.matrix( species[match(rownames(thedata)[species.nodes],species$species),intersect(names(species),colnames(thedata))] )
    thedata[sample.nodes,intersect(names(whales),colnames(thedata))] <- as.matrix( whales[match(rownames(thedata)[sample.nodes],whales$specimen),intersect(names(whales),colnames(thedata))] )
    thedata <- log(thedata)
    havedata <- !is.na(thedata)
} else {
    load(paste(basedir,"thedata-and-covmatrices.Rdata",sep='/'))
}

nvars <- ncol(thedata)

# from parse-correlated-data.R
## tree stuff:
descendants <- get.descendants(tree)  # [i,j] TRUE means that j is below i
n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
rootnode <- which.max(n.offspring)
# edge.indices[k] is the node associated with edge k
# match(j,edge.indices) is the edge associated with node j
edge.indices <- tree$edge[,2] # associate each edge with the downstream node
tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
internal.edges <- setdiff( 1:Nedge(tree), tip.edges )
species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
internal.nodes <- Ntip(tree) + (1:Nnode(tree))

if (!do.simple) {
###### from correlated-traits-analysis.R
# set up with some initial parameters
# the variables are, in order: Length, Testes, Rib-left, Rib-right, Pelvis-left, Pelvis-right

    species.paramnames <- matrix(
                  c('sigmaL','','','','','',
                    'sigmaLdeltaT','betaT','','','','',
                    'sigmaLdeltaP','','sigmaR','','','',
                    'sigmaLdeltaP','','sigmaR','','','',
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


} else {

    ###
    # SIMPLER EXAMPLE

    species.params <- matrix(
                  c('a1','','','','','',
                    'a2','a5','','','','',
                    'a3','','a6','','','',
                    'a3','','a6','','','',
                    'a4','a5','','a7','','',
                    'a4','a5','','a7','',''), nrow=6, byrow=TRUE )
    sample.params <- matrix( 
                  c('b1','','','','','',
                    'b2','b2','','','','',
                    'b3','','b3','-b4','','',
                    'b3','','b3','b4','','',
                    'b4','','','','b5','-b6',
                    'b4','','','','b5','b6'), nrow=6, byrow=TRUE )
    species.varnames <- setdiff( unique( species.params ), '' )
    species.Tind <- which( species.params != '' )  # nonzero elements of species.transmat
    species.Sind <- match( species.params[species.Tind], species.varnames )     # put params[species.Sind] into species.transmat[species.Tind]
    stopifnot( length(species.Sind)==length(species.Tind) )
    sample.justparams <- gsub("^[-+]","",sample.params)
    sample.signs <- ifelse( grepl( "^-", sample.params ), -1, +1 )
    sample.varnames <- setdiff( unique( sample.justparams ), '' )
    sample.Tind <- which( sample.justparams != '' )
    sample.Sind <- match( sample.justparams[sample.Tind], sample.varnames )     # put params[sample.Sind]*sample.Pcoef into sample.transmat[sample.Tind]
    sample.Pcoef <- sample.signs[ sample.Tind ]
    stopifnot( length(sample.Sind)==length(sample.Tind) )

    initpar <- rep(1,max(sample.Sind)+max(species.Sind))
    names(initpar) <- c( paste('a',1:max(species.Sind),sep=''), paste('b',1:max(sample.Sind),sep='') )

    species.params <- initpar[1:max(species.Sind)]
    sample.params <- initpar[max(species.Sind)+1:max(sample.Sind)]
    species.transmat[species.Tind] <- species.params[species.Sind]
    sample.transmat[sample.Tind] <- sample.params[sample.Sind]*sample.Pcoef
    species.covmat <- tcrossprod(species.transmat)
    sample.covmat <- tcrossprod(sample.transmat)
    sptransmat <- species.transmat
    samtransmat <- sample.transmat
### DONE SIMPLER EXAMPLE

}

# [k,j] TRUE if node j is below edge k
edge.descendants <- descendants[edge.indices,]

# normalizing whatnot
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
phylomeans <- lapply( 1:ncol(thedata), function(k) phylomean(thedata[1:Ntip(tree),k], tree=adjtree) )
names(phylomeans) <- colnames(thedata)
tipweights <- lapply( phylomeans, attr, "weights" )
# construct (I-W) term that multiplies the normalized covariance matrix
#  note: indexed by ( variables x tips+nodes )
weightmat <- do.call( cbind, lapply( seq_along(tipweights), function (k) {
        c( rep(0,(k-1)*nrow(thedata)), c(tipweights[[k]],rep(0,Nnode(tree))), rep(0,(ncol(thedata)-k)*nrow(thedata)) )
    } ) )
# here are the W and I-W matrices
weightmat <- weightmat[ , rep(1:ncol(thedata),each=nrow(thedata)) ]
norm.factor <- ( diag( length(thedata) ) - t(weightmat) )
# for likelihood computation, will project data onto the smaller-dimension space:
#  where we have data, and phylomean-centered
center.matrix <- norm.factor[havedata,]
projmatrix.qr <- qr( t(center.matrix) )
projmatrix <- qr.Q( projmatrix.qr )[,1:projmatrix.qr$rank]


###
fake.noise <- function (centered=FALSE) {
    edgediffs <- rnorm( nvars * Nedge(tree) )
    dim(edgediffs) <- c(nvars,Nedge(tree))
    edgediffs <- sweep( edgediffs, 2, sqrt(tree$edge.length), "*" )
    edgediffs[,internal.edges] <- sptransmat %*% edgediffs[,internal.edges] 
    edgediffs[,tip.edges] <- samtransmat %*% edgediffs[,tip.edges] 
    return(edgediffs)
}

###
fake.data <- function (centered=FALSE) {
    edgediffs <- fake.noise(centered=centered)
    # sum down the tree
    simdata <- t( edgediffs %*% edge.descendants )
    if (centered) { 
        simdata <- norm.factor %*% as.vector(simdata)
    } else {
        dimnames(simdata) <- dimnames(thedata)
    }
    return(simdata)
}

simdata <- fake.data()

write.csv(simdata,file='full-simdata.csv',row.names=FALSE)

subsimdata <- simdata
subsimdata[is.na(thedata)] <- NA

write.csv(subsimdata,file='observed-simdata.csv',row.names=FALSE)

# normalize both:
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
sim.phylomeans <- lapply( 1:ncol(subsimdata), function(k) phylomean(subsimdata[1:Ntip(tree),k], tree=adjtree) )
names(sim.phylomeans) <- colnames(subsimdata)
sim.tipweights <- lapply( sim.phylomeans, attr, "weights" )

centered.thedata <- sweep( thedata, 2, sapply(phylomeans, "[[", 1), "-" )
centered.simdata <- sweep( simdata, 2, sapply(sim.phylomeans, "[[", 1), "-" )


###
# look at summary stats

if (FALSE) {

    ##
    # check that covariance matrix matches simulated
    # "internal" edges --
    internal.lengths <- tree$edge.length
    internal.lengths[ tip.edges ] <- 0
    # lengths of shared internal branches leading to the root for all pairs of nodes in the tree
    species.treemat <- shared.branchlengths( tree, internal.lengths, descendants )
    rownames( species.treemat ) <- colnames( species.treemat ) <- c( tree$tip.label, paste("node",Ntip(tree)+1:Nnode(tree),sep='.') )
    dzeros <- diag(species.treemat)==0
    stopifnot( any( species.treemat[dzeros,dzeros] < 1e-8 ) )
    stopifnot(all(abs(cov2cor(species.treemat[!dzeros,!dzeros]))<=1+1e-8))

    # and, in the tips
    tip.lengths <- tree$edge.length
    tip.lengths[ - tip.edges ] <- 0
    sample.treemat <- shared.branchlengths( tree, tip.lengths, descendants )
    rownames( sample.treemat ) <- colnames( sample.treemat ) <- c( tree$tip.label, paste("node",Ntip(tree)+1:Nnode(tree),sep='.') )
    dzeros <- diag(sample.treemat)==0
    stopifnot( any( sample.treemat[dzeros,dzeros] < 1e-8 ) )
    stopifnot(all(abs(cov2cor(sample.treemat[!dzeros,!dzeros]))<=1+1e-8))


    if (!do.simple) {

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

    } else {

        ### SIMPLER EXAMPLE

        make.fullmat <- function (par) {
            # return full covariance matrix for all data (observed and unobserved)
            #  parameters are: ( sigmaL, betaT, betaP, sigmaR, sigmaP ), (zetaL, zetaR, omegaR, zetaP, omegaP), (deltaT, deltaP, deltaR)
            species.params <- par[1:max(species.Sind)]
            sample.params <- par[max(species.Sind)+1:max(sample.Sind)]
            species.transmat[species.Tind] <- species.params[species.Sind] 
            sample.transmat[sample.Tind] <- sample.params[sample.Sind]*sample.Pcoef
            species.covmat <- tcrossprod(species.transmat)
            sample.covmat <- tcrossprod(sample.transmat)
            fullmat <-  kronecker( species.covmat, species.treemat ) + kronecker( sample.covmat, sample.treemat )
            # will want to use 
            # submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
            return( fullmat )
        }

    }

    # projection matrix for omitting all left ribs and pelvises in 'species' and internal nodes
    simple.projmatrix <- diag( length(thedata) )
    simple.havedata <- ! as.vector( (row(thedata) %in% internal.nodes) | ( ( row(thedata) %in% c(species.nodes) ) & ( col(thedata) %in% c(3,5) ) ) )
    simple.projmatrix <- simple.projmatrix[ , simple.havedata ]

    ### END SIMPLER EXAMPLE

    fullmat <- make.fullmat(initpar)
    colnames( fullmat ) <- rownames( fullmat ) <- outer( rownames(thedata), colnames(thedata), paste, sep='.' )

    many.data <- replicate(10000,fake.data())
    dim(many.data) <- c( prod(dim(many.data)[1:2]), dim(many.data)[3] )
    many.cov <- cov(t(many.data))

    ## looks good
    layout(1:2)
    plot( as.vector(fullmat), as.vector(many.cov) )
    abline(0,1)
    plot( as.vector(fullmat), as.vector(fullmat) - as.vector(many.cov) )
    abline(h=0)

    # look more closely
    if (FALSE) {
        many.covs <- list( cov(t(many.data[,1:2500])), cov(t(many.data[,2500+1:2500])), cov(t(many.data[,5000+1:2500])), cov(t(many.data[,7500+1:2500])) )
        plot( as.vector(fullmat), as.vector(fullmat) - as.vector(many.cov), ylim=c(-.2,.2) )
        lapply( seq_along(many.covs), function (k) points( as.vector(fullmat), as.vector(fullmat) - as.vector(many.covs[[k]]), col=rainbow(10)[k] ) )
        points( as.vector(fullmat), as.vector(fullmat) - as.vector(many.cov), ylim=c(-.2,.2), pch=20 )
        abline(h=0)
        plot( as.vector(cov2cor(fullmat)), as.vector(cov2cor(fullmat)) - as.vector(cov2cor(many.cov)), ylim=c(-.05,.05) )
        lapply( seq_along(many.covs), function (k) points( as.vector(cov2cor(fullmat)), as.vector(cov2cor(fullmat)) - as.vector(cov2cor(many.covs[[k]])), col=rainbow(10)[k] ) )
        points( as.vector(cov2cor(fullmat)), as.vector(cov2cor(fullmat)) - as.vector(cov2cor(many.cov)), ylim=c(-.05,.05), pch=20 )
        abline(h=0)
    }

    # projected data?
    proj.data <- crossprod( projmatrix, many.data )
    proj.cov <- cov(t(proj.data))
    proj.fullmat <- crossprod( projmatrix, fullmat %*% projmatrix )
    proj.fchol <- chol(proj.fullmat)
    
    ## again, looks good
    layout(1:2)
    plot( as.vector(proj.fullmat), as.vector(proj.cov) )
    abline(0,1)
    plot( as.vector(cov2cor(proj.fullmat)), as.vector(cov2cor(proj.cov)) )
    abline(0,1)

    # chol transforms to N(0,I)?
    untf.data <- backsolve( proj.fchol, proj.data, transpose=TRUE )
    untf.cov <- cov(t(untf.data))
    range(diag(untf.cov))
    range(untf.cov - diag(diag(untf.cov)))
    hist(untf.data)

    ###
    # check noise?
    many.data <- replicate(10000,fake.data())
    many.noise <- replicate(10000,fake.noise())
    noise.data <- apply( many.noise, 3, function (x) t( x %*%  edge.descendants) )
    dim(noise.data) <- dim(many.data) <- c( prod(dim(many.data)[1:2]), dim(many.data)[3] )
    dim(many.noise) <- c( prod(dim(many.noise)[1:2]), dim(many.noise)[3] )
    noise.cov <- cov(t(many.noise))
    many.cov <- cov(t(many.data))
    noise.data.cov <- cov( t(noise.data) )
    # data covariance agrees?
    range(many.cov-noise.data.cov)
    plot(as.vector(noise.data.cov),as.vector(many.cov)); abline(0,1)
    # noise covariance is as it should be?
    tmp.cov <- matrix(0,nrow=nvars*Nedge(tree),ncol=nvars*Nedge(tree))
    for (k in tip.edges) { tmp.cov[(k-1)*nvars+1:nvars,(k-1)*nvars+1:nvars] <- tree$edge.length[k]*sample.covmat }
    for (k in internal.edges) { tmp.cov[(k-1)*nvars+1:nvars,(k-1)*nvars+1:nvars] <- tree$edge.length[k]*species.covmat }
    # looks good
    range(tmp.cov-noise.cov)
    plot(as.vector(tmp.cov),as.vector(noise.cov))



    ###
    # can use many observations to infer parameters?
    #
    # havedata <- !is.na(thedata)
    havedata <- simple.havedata
    ## unprojected:
    datavec <- many.data[havedata,]  # lots of sim, unprojected
    llfun <- function (par) {
        fullmat <- make.fullmat( par )[havedata,havedata]
        fchol <- chol(fullmat)
        return( sum( backsolve( fchol, datavec, transpose=TRUE )^2 )/2 + ncol(datavec) * sum(log(diag(fchol))) ) 
    }

    ## projected:
    havedata <- !is.na(thedata)
    # datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data
    datavec <- crossprod( projmatrix[havedata,], many.data[havedata,] )  # lots of sims
    ## one simulation
    # simdata <- fake.data(center=TRUE)         # one sim
    # simdata[!havedata] <- NA
    # datavec <- crossprod( projmatrix[havedata,], simdata[havedata] )
    # return negative log-likelihood for gaussian:
    #  parameters are: sigmal, betat, betap, sigmar, sigmap, zetal, zetar, omegar, zetap, omegap, delta
    llfun <- function (par) {
        fullmat <- make.fullmat( par )
        submat <- ( ( crossprod( projmatrix, fullmat) %*% projmatrix ) )
        fchol <- chol(submat)
        return( sum( backsolve( fchol, datavec, transpose=TRUE )^2 )/2 + ncol(datavec) * sum(log(diag(fchol))) ) 
    }
    stopifnot( is.finite(llfun(initpar)) )

    mlestim1 <- optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( maxit=30, fnscale=abs(llfun(initpar)/10) ) )
    # mlestim1 <- optim( par=initpar+runif(length(initpar))/8, fn=llfun, method="Nelder-Mead", control=list( maxit=3000, fnscale=abs(llfun(initpar)/10) ) )
    cbind( rbind(initpar,mlestim1$par), ll=c(llfun(initpar),llfun(mlestim1$par)) )

    # profiles
    xx <- seq(.8, 1.2, length.out=10 )
    layout( matrix(1:16,nrow=4) )
    for (k in 1:length(initpar) ) {
        plot( xx, sapply( xx, function (x) llfun( c( rep(1,k-1), x, rep(1,length(initpar)-k) ) ) ), ylab=names(initpar)[k] )
    }
    # b2-b5 ?
    plot( xx, sapply( xx, function (x) llfun( c( rep(1,8), x, rep(1,2), 2-x, rep(1,1) ) ) ), ylab='b2-b5' )

    # truth
    this.fullmat <- make.fullmat( initpar )
    this.submat <- ( ( crossprod( projmatrix, this.fullmat) %*% projmatrix ) )
    this.fchol <- chol(this.submat)
    this.siginv <- ginv(this.submat)
    # check math:
    all.equal( apply( datavec, 2, function (x) crossprod(x, this.siginv %*% x) ), colSums( backsolve( this.fchol, datavec, transpose=TRUE )^2 ) ) 
    # estimated
    that.fullmat <- make.fullmat( mlestim1$par )
    that.submat <- ( ( crossprod( projmatrix, that.fullmat) %*% projmatrix ) )
    that.fchol <- chol(that.submat)
    # compare to empirical covmatrix
    emp.submat <- cov(t(datavec))
    range(emp.submat-this.submat)
    range(emp.submat-that.submat)
    layout(matrix(1:4,nrow=2))
    plot(as.vector(cov(t(datavec))),this.submat); abline(0,1)
    plot(as.vector(this.submat),as.vector(cov(t(datavec)))-as.vector(this.submat)); abline(h=0)
    plot(as.vector(cov(t(datavec))),that.submat); abline(0,1)
    plot(as.vector(that.submat),as.vector(cov(t(datavec)))-as.vector(that.submat)); abline(h=0)

    # contributions to likelihood?
    this.untf <-  backsolve( this.fchol, datavec, transpose=TRUE )
    that.untf <-  backsolve( that.fchol, datavec, transpose=TRUE )
    summary(as.vector(this.untf))
    summary(as.vector(that.untf))
    hist( this.untf )
    qqnorm(this.untf); abline(0,1)
    hist( that.untf ) 
    qqnorm(that.untf); abline(0,1)
    c( this=sum( this.untf^2 ), that=sum(that.untf^2) )
    ncol(datavec) * c( this=sum(log(diag(this.fchol))), that=sum(log(diag(that.fchol))) )


    ####
    # fit by matching covariance matrix?

    pmat <- projmatrix
    # pmat <- simple.projmatrix
    # obs.cov <- crossprod( pmat, many.cov ) %*% pmat
    obs.cov <- crossprod( pmat, many.cov ) %*% pmat
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )

    covfit <- function (par) {
        # fullmat <- make.fullmat( par )[havedata,havedata]
        # fchol <- chol(fullmat)
        fullmat <- make.fullmat( par )
        # submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
        submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
        fchol <- chol(submat)
        return( sum( ( obs.cov - submat )^2 ) )
    }
    stopifnot( is.finite(covfit(initpar)) )

    moment.estim1 <- optim( par=initpar, fn=covfit, method="Nelder-Mead", control=list( maxit=30, fnscale=abs(covfit(initpar)/10) ) )
    # moment.estim1 <- optim( par=initpar+runif(length(initpar))/8, fn=covfit, method="Nelder-Mead", control=list( maxit=3000 ) )
    cbind( rbind(initpar,moment.estim1$par), ll=c(llfun(initpar),llfun(moment.estim1$par)), mom=c(covfit(initpar),covfit(moment.estim1$par)) )

    est.cov <- make.fullmat( moment.estim1$par )
    est.subcov <- ( ( crossprod( pmat, est.cov) %*% pmat ) )
    internal.cols <- 1 + ( row(thedata) %in% internal.nodes ) + ( col(thedata) %in% internal.nodes )
    layout(matrix(1:4,nrow=2))
    plot( as.vector(est.cov), as.vector(many.cov), xlab='fitted values', ylab='observed', col=col(thedata) ); abline(0,1)
    plot( as.vector(fullmat), as.vector(many.cov), xlab='true values', ylab='observed', col=col(thedata) ); abline(0,1)
    plot( as.vector(est.subcov), as.vector(obs.cov), xlab='fitted values', ylab='observed', col=col(thedata) ); abline(0,1)
    plot( as.vector(submat), as.vector(obs.cov), xlab='true values', ylab='observed', col=col(thedata) ); abline(0,1)

    # profiles
    xx <- seq(.8, 1.2, length.out=25 )
    layout( matrix(1:16,nrow=4) )
    for (k in 1:length(initpar) ) {
        plot( xx, sapply( xx, function (x) covfit( c( rep(1,k-1), x, c(rep(1,length(initpar)-k)) ) ) ) )
    }


    ######
    ## look at predicted sum-square-differences between paired bones

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

    true.bonediffs <- bonediffmat %*% ifelse( is.na(as.vector(thedata)), 0, thedata )
    sim.bonediffs <- bonediffmat %*% ifelse( is.na(as.vector(thedata)), 0, subsimdata )
    mean(true.bonediffs); sqrt( sum(true.bonediffs^2) )
    mean(sim.bonediffs); sqrt( sum(sim.bonediffs^2) )

    ## compare differences across phylogenetic distances
    distmat <- treedist( tree, edge.length=ifelse(seq_along(tree$edge.length)%in%tip.edges,0,tree$edge.length) )
    dimnames(distmat) <- list( rownames(thedata), rownames(thedata) )
    layout(matrix(1:12,nrow=2))
    for (k in 1:ncol(thedata)) { plot( as.vector( distmat ), outer(thedata[,k],thedata[,k],"-"), main=colnames(thedata)[k] ) }
    for (k in 1:ncol(simdata)) { plot( as.vector( distmat ), outer(simdata[,k],simdata[,k],"-"), main=colnames(simdata)[k] ) }



    ## predicted sum-square-differences in closely related clades


}

