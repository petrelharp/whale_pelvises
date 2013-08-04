#!/usr/bin/R

# implements things discussed in correlated-trats.tex

require(ape)

phylomean <- function (x,tree) {
    # "phylogenetic" mean of values given at the tips.
    # this must be implemented in ape somewhere.
    #  as in Felsenstein '73.
    # return also the vector of weights so that the answer is sum(phylomean*x)
    if (!is.null(names(x))) {
        tree <- drop.tip( tree, setdiff( tree$tip.label, names(x) ) )
        if (!setequal(names(x),tree$tip.label)) {
            warning("Names of x and tip labels don't match.")
        }
        xperm <- match( names(x), tree$tip.label )
        x <- x[ xperm ]
    }
    tree <- drop.tip( tree, tree$tip.label[is.na(x)] )
    meanvals <- c(x[!is.na(x)],rep(NA,Nnode(tree)))
    thetips <- 1:Ntip(tree)
    weights <- matrix(0, nrow=Nnode(tree)+Ntip(tree), ncol=length(x))
    colnames(weights) <- names(x)
    for (k in 1:Ntip(tree)) { weights[k,match(k,cumsum(!is.na(x)))] <- 1 }
    theroot <- setdiff(tree$edge[,1],tree$edge[,2])
    edgelens <- tree$edge.length
    for (k in 1:Nnode(tree)) {
        tip.offspring <- tapply( tree$edge[,2], tree$edge[,1], function (y) all(y %in% thetips) )
        cherry <- as.numeric( names(tip.offspring)[tip.offspring] )[1] # node we will prune to
        stems <- which(tree$edge[,1]==cherry)  # edges to be pruned
        chstem <- which(tree$edge[,2]==cherry) # edge leading to basal node
        chtips <- tree$edge[stems,2]  # nodes at ends of to-be-pruned edges
        chweights <- (1/edgelens[stems]) / sum(1/edgelens[stems])
        meanvals[cherry] <- sum( meanvals[chtips] * chweights )
        weights[cherry,] <- colSums( weights[chtips,] * chweights )
        edgelens[chstem] <- edgelens[ chstem ] + 1/sum( 1 / edgelens[stems] )
        thetips <- c( cherry, setdiff( thetips, chtips ) )
        # cbind( meanvals, rowSums( weights * x[col(weights)], na.rm=TRUE ) )
    }
    stopifnot(cherry==theroot)
    ans <- meanvals[theroot] 
    stopifnot( all.equal( ans, sum(weights[theroot,]*x,na.rm=TRUE), check.attributes=FALSE ) )
    attr(ans,"weights") <- weights[theroot,order(xperm)]
    return( ans )
}

get.descendants  <- function (tree) {
    # return the node-by-node matrix with [i,j] TRUE meaning that j is below i
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
    return(descendants)
}

tip.to.edge.weights <- function (w,tree,descendants=get.descendants(tree)) {
    # given a set of weights w on the tips of the rooted tree,
    # return edge weights: the weight on an edge is the sum of all weights on descendant tips.
    if (!is.null(names(w)) & all( names(w) %in% tree$tip.label ) ) { w <- w[tree$tip.label] }
    if (length(w) != Ntip(tree)) { stop("w should be named or have as many tips as tree") }
    edge.indices <- tree$edge[,2] # associate each edge with the downstream node
    descendants <- descendants[ edge.indices, 1:Ntip(tree) ]
    return( rowSums( w[col(descendants)] * descendants, na.rm=TRUE ) )
}

shared.branchlengths <- function (tree,edge.length=tree$edge.length,descendants=get.descendants(tree)) {
    # pairwise lengths of shared internal branches leading to the root for all pairs of nodes in the tree
    edge.indices <- tree$edge[,2] # associate each edge with the downstream node
    tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
    n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
    rootnode <- which.max(n.offspring)
    treedist <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) )
    for (j in 1:nrow(treedist))  for (k in 1:ncol(treedist)) {
        jdesc <- descendants[,j] 
        kdesc <- descendants[,k] 
        contributes <- match( which( ( jdesc & kdesc ) ), c(edge.indices,rootnode) )
        treedist[j,k] <- sum( c(edge.length,0)[ contributes ] )
    }
    return( treedist )
}

treedist <- function (tree,edge.length=tree$edge.length,descendants=get.descendants(tree)) {
    # pairwise lengths of shared internal branches leading to the root for all pairs of nodes in the tree
    edge.indices <- tree$edge[,2] # associate each edge with the downstream node
    tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
    n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
    rootnode <- which.max(n.offspring)
    treedist <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) )
    for (j in 1:nrow(treedist))  for (k in 1:ncol(treedist)) {
        jdesc <- descendants[,j] 
        kdesc <- descendants[,k] 
        contributes <- match( which( ( jdesc | kdesc ) & ! ( jdesc & kdesc ) ), c(edge.indices,rootnode) )
        treedist[j,k] <- sum( c(edge.length,0)[ contributes ] )
    }
    return( treedist )
}

predgaus <- function( obsvec, meanvec, covmat ) {
    # fill in NAs in 'obsvec' by kriging
    ii <- is.na(obsvec)
    krig <- meanvec[ ii ] + covmat[ ii, !ii ] %*% solve( covmat[!ii,!ii], (obsvec-meanvec)[!ii] ) 
    obsvec[ii] <- krig
    return(obsvec)
}
