#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(Matrix)

runfiles <- list.files(path="shape-likelihood-surface",pattern=".*RData",full.names=TRUE)
runs <- lapply( runfiles, function (x) { tmpenv <- environment(); load(x,envir=tmpenv); as.list(tmpenv) } )

pargrid <- runs[[1]]$pargrid
stopifnot( all( sapply( runs, function (x) all.equal(pargrid, x$pargrid) ) ) )
sigma2S.vals <- sort( unique( pargrid$sigma2S ) )
gammaP.vals <- sort( unique( pargrid$gammaP ) )

pelvic.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "pelvic.sgrids" ) ) )
sub.pelvic.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "sub.pelvic.sgrids" ) ) )
rib.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "rib.sgrids" ) ) )

## normalize, within each
pelvic.normed <- exp(pelvic.sgrids-max(pelvic.sgrids)) 
pelvic.normed <- sweep( pelvic.normed, 2, colSums(pelvic.normed), "/" )
rib.normed <- exp(rib.sgrids-max(rib.sgrids)) 
rib.normed <- sweep( rib.normed, 2, colSums(rib.normed), "/" )
sub.pelvic.normed <- exp(sub.pelvic.sgrids-max(sub.pelvic.sgrids)) 
sub.pelvic.normed <- sweep( sub.pelvic.normed, 2, colSums(sub.pelvic.normed), "/" )

pargrid$pelvic <- rowMeans( pelvic.normed )
pargrid$rib <- rowMeans( rib.normed )
pargrid$sub.pelvic <- rowMeans( sub.pelvic.normed )
scalefac <- 500

##
# marginal likelihoods of ks
sapply( c("pelvic","sub.pelvic","rib"), function (x) tapply( pargrid[[x]]/sum(pargrid[[x]]), pargrid$ks, sum ) )

# joint distributions:
pdf(file="joint-shape-posteriors.pdf",width=8,height=8,pointsize=10)
layout(matrix(1:4,nrow=2))
for (this.ks in sort(unique(pargrid$ks))) {
    with( subset( pargrid, ks==this.ks ), {
                plot( sigma2S, gammaP, cex=scalefac*pelvic, pch=21, col='grey', bg=adjustcolor('red',.3), main=paste("ks = ", this.ks) )
                points( sigma2S, gammaP, cex=scalefac*rib, pch=21, col='grey', bg=adjustcolor('blue',.3) )
                points( sigma2S, gammaP, cex=scalefac*sub.pelvic, pch=21, col='grey', bg=adjustcolor('grey',.3) )
            } )
    if (this.ks==min(pargrid$ks)) { legend("topleft",pch=21,pt.bg=adjustcolor(c('red','blue','grey'),.3), pt.cex=2, legend=c("pelvic","rib","pelvic, complete obs") ) }
}
dev.off()

pdf(file="hist-shape-posteriors.pdf",width=7,height=5,pointsize=10)
layout(t(1:2))
gP.dists <- lapply( c("pelvic","sub.pelvic","rib"), function (type) {
            tapply( pargrid[[type]], pargrid$gammaP, sum ) } )
plot( 0, type='n', xlab='gammaP', ylab='probability', xlim=range(pargrid$gammaP), ylim=range(gP.dists) )
fillcols <- adjustcolor(c('red','blue','grey'),.3)
for (k in 1:3) {
    polygon( c( min(gammaP.vals), gammaP.vals, max(gammaP.vals) ), c(0,gP.dists[[k]],0), col=fillcols[k] )
}
sS.dists <- lapply( c("pelvic","sub.pelvic","rib"), function (type) {
            tapply( pargrid[[type]], pargrid$sigma2S, sum ) } )
plot( 0, type='n', xlab='sigma2S', ylab='probability', xlim=range(pargrid$sigma2S), ylim=range(sS.dists) )
fillcols <- adjustcolor(c('red','blue','grey'),.3)
for (k in 1:3) {
    polygon( c( min(sigma2S.vals), sigma2S.vals, max(sigma2S.vals) ), c(0,sS.dists[[k]],0), col=fillcols[k] )
}
dev.off()

if (interactive()) {
    # look at variability across inferred testes size
    layout(1:3)
    for (type in c("pelvic","sub.pelvic","rib")) {
        this.sgrids <- get( paste(type,".sgrids",sep='') )
        usethese <- (pargrid$ks==20)
        sigma2S.vals <- sort( unique( pargrid$sigma2S ) )
        gammaP.vals <- sort( unique( pargrid$gammaP ) )
        z <- exp(this.sgrids[,1][usethese])
        dim(z) <- c( length(sigma2S.vals), length(gammaP.vals) )
        contour( sigma2S.vals, gammaP.vals, z, col=NA, main=type )
        for (this.ks in sort(unique(pargrid$ks))) {
            usethese <- (pargrid$ks==this.ks)
            for (k in sample(ncol(this.sgrids),30)) {
                z <- exp(this.sgrids[,k][usethese])
                dim(z) <- c( length(sigma2S.vals), length(gammaP.vals) )
                contour( sigma2S.vals, gammaP.vals, z, add=TRUE, col=adjustcolor(this.ks-16,.2) )
            }
        }
    }
}
