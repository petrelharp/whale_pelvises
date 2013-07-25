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
# associate P and Q with each internal branch of the tree:
#  these parameters are: theta = sigmaL, betaT, betaP, sigmaR, sigmaP.
#  The matrix P is such that P %*% theta gives the NONZERO elements of the corresponding sqrt-covariance matrix.
#  ... which are (branch length) * (sigmaL, sigmaL, sigmaL, sigmaL, betaT, betaP, sigmaR, sigmaP)
#  and the matrix Q is similar, except says where to put delta (on all sigmaL but the first one)
species.Pmat <- cbind( c(1,1,1,1,0,0,0,0),
           c(0,0,0,0,1,0,0,0),
           c(0,0,0,0,0,1,0,0),
           c(0,0,0,0,0,0,1,0),
           c(0,0,0,0,0,0,0,1) )
species.Qmat <- c(0,1,1,1,0,0,0,0)

# Now add tips for samples:
# When constructing P, Q, 
#  the paramters are: zetaL, zetaR, omegaR, zetaP, omegaP
# and the nonzero elements, in order, are
#  (zetaL, zetaL, zetaL, zetaL, zetaL, zetaR, zetaR, omegaR, -omegaR, zetaP, zetaP, omegaP, -omegaP)
# again, put delta on all zetaL but the first one
#  ... same for every sample edge
sample.Pmat <- cbind( c(1,1,1,1,1,0,0,0,0,0,0,0,0),
                            c(0,0,0,0,0,1,1,0,0,0,0,0,0),
                            c(0,0,0,0,0,0,0,1,-1,0,0,0,0),
                            c(0,0,0,0,0,0,0,0,0,1,1,0,0),
                            c(0,0,0,0,0,0,0,0,0,0,0,1,-1) )
sample.Qmat <- c(0,1,1,1,1,0,0,0,0,0,0,0,0)
for (sp in intersect(levels(whales$species),tree$tip.label)) {
    theseones <- (whales$species==sp)
    nsamps <- sum(theseones)
    tmptree <- list(
        edge = cbind( rep(nsamps+1,nsamps), 1:nsamps ),
        edge.length = rep(within_length,nsamps),
        Nnode = 1,
        tip.label = levels(whales$specimen)[as.numeric(whales$specimen[theseones])],
        root.edge = within_length
        )
    class(tmptree) <- "phylo"
    tree <- bind.tree(x=tree,y=tmptree,where=match(sp,tree$tip.label))
}

tree <- drop.tip( tree, setdiff( tree$tip.label, whales$specimen ) )
stopifnot( setequal( whales$specimen, tree$tip.label ) )

# construct full matrix
nvars <- 4
full.covmat <- matrix(0,nrow=nvars*Ntip(tree),ncol=nvars*Ntip(tree))
rowtip <- ( ( row(full.covmat)-1) %/% nvars ) + 1 
coltip <- ( ( col(full.covmat)-1) %/% nvars ) + 1 

# pre-compute the number of times each edge's covariance matrix enters into the full covariance matrix
# as well as the vector indicating where it goes
adjacency <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) ) 
adjacency[ tree$edge ] <- 1  # adjacency matrix: [i,j] means that node j is directly below node i
descendants <- apower <- adjacency  # [i,j] means that node j is somewhere below node i
while ( any(apower>0) ) {
    apower <- apower %*% adjacency
    descendants <- descendants + apower
}
n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
rootnode <- which.max(n.offspring)
edge.indices <- setdiff( 1:Nnode(tree), rootnode ) # associate each edge with the downstream node

# set up with some initial parameters
species.params <- c( sigmaL=.2, betaT=.1, betaP=.3, sigmaR=.25, sigmaP=.05 )
sample.params <- c( zetaL=.1, zetaR=.1, omegaR=.1, zetaP=.1, omegaP=.2 )
delta <- 1.2
species.transmat <- Matrix( c(1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,1), nrow=4, sparse=TRUE )
sample.transmat <- Matrix( c(1,1,1,1,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,1,1), nrow=5, sparse=TRUE )
species.transmat@x <- as.vector( ( species.Pmat %*% species.params ) * ( species.Qmat * delta ) )
sample.transmat@x <- as.vector( ( sample.Pmat %*% sample.params ) * ( sample.Qmat * delta ) )

# list corresponding to each edge
covmats <- c( 
        lapply( species_tree$edge.length, function (len) { len * tcrossprod( species.transmat ) } ),
        lapply( 1:(Nedge(tree)-Nedge(species_tree)), function (k) { tcrossprod( sample.transmat ) } )
    )

for (k in 1:Nedge(tree)) {
    # add covmat[[k]] to the appropriate pieces of full.covmat

}


##

