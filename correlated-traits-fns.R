#!/usr/bin/R

# implements things discussed in correlated-trats.tex

require(ape)

phylomean <- function (x,tree) {
    # this must be implemented in ape somewhere.
    #  as in Felsenstein '73.
    # return also the vector of weights so that the answer is sum(phylomean*x)
    if (!is.null(names(x))) {
        tree <- drop.tip( tree, setdiff( tree$tip.label, names(x) ) )
        if (!setequal(names(x),tree$tip.label)) {
            warning("Names of x and tip labels don't match.")
        }
        x <- x[tree$tip.label]
    }
    tree <- drop.tip( tree, tree$tip.label[is.na(x)] )
    meanvals <- c(x[!is.na(x)],rep(NA,Nnode(tree)))
    thetips <- 1:Ntip(tree)
    weights <- matrix(0, nrow=Nnode(tree)+Ntip(tree), ncol=length(x))
    for (k in thetips) { weights[k,cumsum(!is.na(x))[k]] <- 1 }
    theroot <- setdiff(tree$edge[,1],tree$edge[,2])
    edgelens <- tree$edge.length
    for (k in 1:Nnode(tree)) {
        tip.offspring <- tapply( tree$edge[,2], tree$edge[,1], function (y) sum(y %in% thetips) )
        cherry <- as.numeric( names(tip.offspring)[tip.offspring==2] )[1]
        stems <- which(tree$edge[,1]==cherry)
        chstem <- which(tree$edge[,2]==cherry)
        chtips <- tree$edge[stems,2]
        chweights <- (1/edgelens[stems]) / sum(1/edgelens[stems])
        meanvals[cherry] <- sum( meanvals[chtips] * chweights )
        weights[cherry,] <- ( weights[chtips[1],] * chweights[1] + weights[chtips[2],] * chweights[2] )
        edgelens[chstem] <- edgelens[ chstem ] + 1/sum( 1 / edgelens[stems] )
        thetips <- c( cherry, setdiff( thetips, tree$edge[stems,2] ) )
    }
    stopifnot(cherry==theroot)
    stopifnot(meanvals[theroot] == sum(weights[theroot,]*x) )
    ans <- meanvals[theroot] 
    attr(ans,"weights") <- weights[theroot,]
    return( ans )
}

