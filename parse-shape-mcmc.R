

# all the data
layout(1:3)
for (k in 1:3) {
    tmp <- hist(c(pelvic.mcrun$batch[,k],rib.mcrun$batch[,k]),breaks=50,plot=FALSE)
    hist(rib.mcrun$batch[,k],breaks=50,col=adjustcolor('blue',.5),main=names(initpar)[k],xlim=range(tmp$breaks))
    hist(pelvic.mcrun$batch[,k],breaks=50,col=adjustcolor('red',.5),add=TRUE)
}

# complete obs
layout(1:3)
for (k in 1:3) {
    tmp <- hist(c(sub.pelvic.mcrun$batch[,k],sub.rib.mcrun$batch[,k]),breaks=50,plot=FALSE)
    hist(sub.rib.mcrun$batch[,k],breaks=50,col=adjustcolor('blue',.5),main=names(initpar)[k],xlim=range(tmp$breaks))
    hist(sub.pelvic.mcrun$batch[,k],breaks=50,col=adjustcolor('red',.5),add=TRUE)
}

if (FALSE) {

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

    # moment estimate
    have.dpelvic <- !is.na(pelvicdiff)
    have.both <- !is.na(pelvicdiff) & !is.na(ribdiff)
    # usethese <- have.both
    usethese <- have.dpelvic
    datavec <- pelvicdiff[usethese] / initpar[1]
    f1 <- function (spar) {
        sigma2S <- spar[1]
        gammaP <- spar[2]
        xi2P <- spar[3]
        elens <- ( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
        treemat <- treedist( adjtree, edge.length=elens, descendants=descendants )
        return( sum( ( datavec - treemat[usethese] )^2 ) )
    }

    est.vals <- optim( par=initpar[2:4], fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=100) )
    est.vals <- optim( par=est.vals$par, fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=10) )
    est.vals <- optim( par=est.vals$par, fn=f1, method="BFGS", control=list(fnscale=1e4,maxit=1000) )
    # Using everything:
    # only took 31 iterations
    #    sigma2S     gammaP       xi2P 
    # 0.11451506 0.02237210 0.03384165 
    #
    # Using only ones with ribs also:
    #    sigma2S     gammaP       xi2P 
    # 0.09869252 0.01409928 0.03361052 

    new.treemat <- treedist( adjtree, edge.length=scale.brlens(est.vals$par), descendants=descendants )

    if (interactive()) {
        layout(t(1:2))
        plot( treemat[usethese], datavec )
        abline(0,1)
        plot( new.treemat[usethese], datavec )
        abline(0,1)
    }
}

if (FALSE) {

##############
# pseudolikelihood NOT WORKING
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

}
