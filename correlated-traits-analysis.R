#!/usr/bin/R

source("correlated-traits-fns.R")
require(Matrix)

tree_file <- "consensusTree_ALL_CETACEA.tree"
species_tree<-read.nexus(file=tree_file)

bones <- read.table("50_make_datamatrix.out", header=TRUE)

# rearrange to have one entry per whale.
whichbone <- factor( paste(bones$side, bones$bone, sep='.') )
usevars <- setdiff(names(bones),c("side","bone"))
by.bone <- tapply( 1:nrow(bones), whichbone, function (k) bones[k,,drop=FALSE] )
for (k in seq_along(by.bone)) {
    by.bone[[k]] <- by.bone[[k]][usevars]
    names(by.bone[[k]])[match("absolute_volume",names(by.bone[[k]]))] <- names(by.bone)[k]
}
whales <- by.bone[[1]]
for (k in 2:length(by.bone)) { whales <- merge( whales, by.bone[[k]], all=TRUE ) }

# set up the tree
within_length <- 1
tree <- species_tree
for (sp in intersect(levels(whales$species),tree$tip.label)) {
    theseones <- (whales$species==sp)
    nsamps <- sum(theseones)
    tmptree <- list(  # include zero-length edge (to first node) for species observation
        edge = cbind( rep(nsamps+2,nsamps+1), 1:(nsamps+1) ),
        edge.length = c(0,rep(within_length,nsamps)),
        Nnode = 1,
        tip.label = c(sp,levels(whales$specimen)[as.numeric(whales$specimen[theseones])]),
        root.edge = 0  # this is not a length.
        )
    class(tmptree) <- "phylo"
    tree <- bind.tree(x=tree,y=tmptree,where=match(sp,tree$tip.label))
}

tree <- drop.tip( tree, setdiff( tree$tip.label, c(levels(whales$specimen),levels(whales$species)) ) )
stopifnot( setequal( c(levels(whales$specimen),levels(whales$species)), tree$tip.label ) )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% 1:Ntip(tree) ] %in% c(0,within_length) ) )

bones <- read.table("50_make_datamatrix.out", header=TRUE)

# pre-compute the number of times each edge's covariance matrix enters into the full covariance matrix
# as well as the vector indicating where it goes
adjacency <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) ) 
adjacency[ tree$edge ] <- 1  # adjacency matrix: [i,j] means that node j is directly below node i
# descendants is indexed by NODES
descendants <- apower <- adjacency  # [i,j] means that node j is somewhere below node i
while ( any(apower>0) ) {
    apower <- apower %*% adjacency
    descendants <- descendants + apower
}
stopifnot( all( descendants %in% c(0,1) ) )
descendants <- ( descendants > 0 )
diag(descendants) <- TRUE
n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
rootnode <- which.max(n.offspring)
# edge.indices[k] is the node associated with edge k
# match(j,edge.indices) is the edge associated with node j
edge.indices <- tree$edge[,2] # associate each edge with the downstream node
tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% species.nodes ] == 0 ) )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% sample.nodes ] == within_length ) )
stopifnot( ! ( rootnode %in% edge.indices ) )
# # construct the (edges) x (tips) x (tips) array that says which branches lie between which tips
# contributes.array <- array( FALSE, dim=c(Nedge(tree), Ntip(tree),Ntip(tree)) )
# for (j in 1:Ntip(tree)) for (k in 1:Ntip(tree)) {
#     jdesc <- descendants[,j] 
#     kdesc <- descendants[,k] 
#     contributes.array[,j,k][ match( which( ( jdesc | kdesc ) & ! ( jdesc & kdesc ) ), edge.indices ) ] <- TRUE
# }
# dimnames( contributes ) <- list( edge.indices, tree$tip.label, tree$tip.label )
# pairwise distances between all nodes in the tree
treedist <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) )
internal.lengths <- tree$edge.length
internal.lengths[ tip.edges ] <- 0
for (j in 1:nrow(treedist))  for (k in 1:ncol(treedist)) {
    jdesc <- descendants[,j] 
    kdesc <- descendants[,k] 
    contributes <- match( which( ( jdesc | kdesc ) & ! ( jdesc & kdesc ) ), c(edge.indices,rootnode) )
    treedist[j,k] <- sum( c(internal.lengths,0)[ contributes ] )
}
# and, "dists" along tips
tipdist <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) )
for (j in 1:nrow(tipdist))  for (k in 1:ncol(tipdist)) {
    if (j!=k) {
        jdesc <- descendants[,j] 
        kdesc <- descendants[,k] 
        tipdist[j,k] <- sum( c(j,k) %in% sample.nodes )  # number of j,k that are sample tips
    }
}
stopifnot( abs( cophenetic.phylo(tree) - (treedist+tipdist)[1:Ntip(tree),1:Ntip(tree)] ) < 1e-8 )

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
species.params <- c( sigmaL=.2, betaT=.1, betaP=.3, sigmaR=.25, sigmaP=.05 )
sample.params <- c( zetaL=.1, zetaR=.1, omegaR=.1, zetaP=.1, omegaP=.2 )
delta <- 1.2
species.transmat <- Matrix( 
              c(1,0,0,0,
                1,1,0,0,
                1,0,1,0,
                1,0,1,0,
                1,1,0,1,
                1,1,0,1), nrow=6, byrow=TRUE, sparse=TRUE )
sample.transmat <- Matrix( 
              c(1,1,1,1,1,
                0,0,0,0,0,
                1,1,1,0,0,
                1,1,-1,0,0,
                1,0,0,1,1,
                1,0,0,1,-1), nrow=6, byrow=TRUE, sparse=TRUE )
species.transmat@x <- as.vector( ( species.params[species.Pmat] ) * ( 1 + species.Qmat * (delta-1) ) )
sample.transmat@x <- as.vector( ( sample.params[sample.Pmat] * sample.Pcoef ) * ( 1 + sample.Qmat * (delta-1) ) )
species.covmat <- tcrossprod(species.transmat)
sample.covmat <- tcrossprod(sample.transmat)

# construct full matrix
fullmat <- kronecker( treedist, species.covmat ) + kronecker( tipdist, sample.covmat )
