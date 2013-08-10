#!/usr/bin/R

require(ape)
require(Matrix)

basedir <- scriptdir <- ".."
source(paste(scriptdir,"correlated-traits-fns.R",sep="/"))

tree_file <- paste(basedir, "consensusTree_ALL_CETACEA.tree", sep='/')
species_tree<-read.nexus(file=tree_file)
load( paste(basedir,"all-sample-tree.RData",sep='/') ) # gets tree

load(paste(basedir,"thedata-and-covmatrices.Rdata",sep='/'))
whales <- read.csv(paste(basedir,"whales.csv",sep='/'),header=TRUE)

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

###### from correlated-traits-analysis.R
# associate P and Q with each internal branch of the tree:
#  these parameters are: theta = sigmaL, betaT, betaP, sigmaR, sigmaP.
#  The vector Pmat is such that theta[Pmat] gives the NONZERO elements of the corresponding sqrt-covariance matrix.
#  ... which are (branch length) * (sigmaL, sigmaL, sigmaL, sigmaL, betaT, betaP, sigmaR, sigmaP)
#  and the matrix Q is similar, except says where to put delta (on all sigmaL but the first one)
species.Pmat <- c(1,1,1,1,1,1,2,3,3,4,4,5,5)
species.Qmat <- c(0,1,1,1,1,1,0,0,0,0,0,0,0)
species.delta.Pmat <- c(0,1,2,2,3,3,0,0,0,0,0,0,0)

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
sample.delta.Pmat <- c(0,2,2,3,3,0,0,0,0,0,0,0,0)

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
        deltaT=sqrt(1.4),
        deltaP=sqrt(1.3),
        deltaR=sqrt(1.2)
    )

species.params <- initpar[1:5]
sample.params <- initpar[5+1:5]
delta <- initpar[11:13]
species.transmat@x <- as.vector( ( species.params[species.Pmat] ) * c(1,delta)[1+species.delta.Pmat] )
sample.transmat@x <- as.vector( ( sample.params[sample.Pmat] * sample.Pcoef ) * c(1,delta)[1+sample.delta.Pmat] )
species.covmat <- as.matrix( tcrossprod(species.transmat) )
sample.covmat <- as.matrix( tcrossprod(sample.transmat) )
sptransmat <- cbind( as.matrix(species.transmat), matrix(0,ncol=2,nrow=6) )
samtransmat <- cbind( as.matrix(sample.transmat), matrix(0,ncol=1,nrow=6) )

###
edgediffs <- rnorm( nvars * Nedge(tree) )
dim(edgediffs) <- c(nvars,Nedge(tree))
edgediffs[,internal.edges] <- sptransmat %*% edgediffs[,internal.edges] 
edgediffs[,tip.edges] <- samtransmat %*% edgediffs[,tip.edges] 

# sum down the tree
# [k,j] TRUE if node j is below edge k
edge.descendants <- descendants[edge.indices,]
simdata <- t( edgediffs %*% edge.descendants )
dimnames(simdata) <- dimnames(thedata)

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

