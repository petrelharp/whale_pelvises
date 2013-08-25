#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)

####
# update branch lengths function
# internal branches setup
internal.lengths <- tree$edge.length
internal.lengths[ tip.edges ] <- 0
# and, the tips
tip.lengths <- tree$edge.length
tip.lengths[ - tip.edges ] <- 0

# given parameters get rescaled branch lengths
scale.brlens <- function (par) {
    # par = sigma2S, gammaP, xi2P
    sigma2S <- par[1]
    gammaP <- par[2]
    xi2P <- par[3]
    return( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
}

treemat <- treedist( adjtree, edge.length=scale.brlens(initpar[2:4]), descendants=descendants )

# moment estimate?
have.dpelvic <- !is.na(pelvicdiff)
datavec <- pelvicdiff[have.dpelvic] / initpar[1]
f1 <- function (spar) {
    sigma2S <- spar[1]
    gammaP <- spar[2]
    xi2P <- spar[3]
    elens <- ( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
    treemat <- treedist( adjtree, edge.length=elens, descendants=descendants )
    return( sum( ( datavec - treemat[have.dpelvic] )^2 ) )
}

est.vals <- optim( par=initpar[2:4], fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=10) )
est.vals <- optim( par=est.vals$par, fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=10) )

new.treemat <- treedist( adjtree, edge.length=scale.brlens(est.vals$par), descendants=descendants )

if (interactive()) {
    layout(t(1:2))
    plot( treemat[have.dpelvic], datavec )
    abline(0,1)
    plot( new.treemat[have.dpelvic], datavec )
    abline(0,1)
}

# pseudolikelihood
have.dpelvic <- !is.na(pelvicdiff)
chisq.datavec <- pelvicdiff[have.dpelvic]
f2 <- function (spar) {
    sigma2S <- spar[1]
    gammaP <- spar[2]
    xi2P <- spar[3]
    elens <- ( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
    treemat <- treedist( adjtree, edge.length=elens, descendants=descendants )
    return( (-1)*sum( dchisq( chisq.datavec/treemat[have.dpelvic], initpar[1], log=TRUE ) ) )
}

est.vals <- optim( par=initpar[2:4], fn=f2, method="BFGS", control=list(fnscale=1e4,maxit=10) )
est.vals <- optim( par=est.vals$par, fn=f2, method="BFGS", control=list(fnscale=1e4,maxit=10) )

########
# add tips
rib.lengths <- tree$edge.length
rib.lengths[ tip.edges ] <- xi2R
pelvic.lengths <- tree$edge.length
pelvic.lengths[ tip.edges ] <- xi2P
rib.treemat <- treedist( tree, edge.length=rib.lengths )
pelvic.treemat <- treedist( tree, edge.length=pelvic.lengths )

layout(matrix(1:4,nrow=2))
plot( as.vector(rib.treemat), as.vector(ribdiff), ylab="rib differences" )
plot( as.vector(pelvic.treemat), as.vector(pelvicdiff), ylab="pelvis differences" )

