#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(Matrix)

runfiles <- list.files(path="shape-likelihood-surface",pattern=".*RData",full.names=TRUE)
runs <- lapply( runfiles, function (x) { tmpenv <- environment(); load(x,envir=tmpenv); as.list(tmpenv) } )

pargrid <- runs[[1]]$pargrid
stopifnot( all( sapply( runs, function (x) all.equal(pargrid, x$pargrid) ) ) )

pelvic.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "pelvic.sgrids" ) ) )
sub.pelvic.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "sub.pelvic.sgrids" ) ) )
rib.sgrids <- do.call( cbind, do.call( c, lapply( runs, "[[", "rib.sgrids" ) ) )

pargrid$pelvic <- rowSums(exp(pelvic.sgrids+120))
pargrid$sub.pelvic <- rowSums(exp(sub.pelvic.sgrids+20))
pargrid$rib <- rowSums(exp(rib.sgrids-40))
# for ( x in c("pelvic","sub.pelvic","rib")) { pargrid[[x]] <- pargrid[[x]] / sum(pargrid[[x]]) }

##
# marginal likelihoods of ks
sapply( c("pelvic","sub.pelvic","rib"), function (x) tapply( pargrid[[x]]/sum(pargrid[[x]]), pargrid$ks, sum ) )

layout(matrix(1:4,nrow=2))
scalefac <- 10
for (this.ks in sort(unique(pargrid$ks))) {
    with( subset( pargrid, ks==this.ks ), {
                plot( sigma2S, gammaP, cex=scalefac*pelvic, pch=21, col='grey', bg=adjustcolor('red',.3), main=paste("ks = ", this.ks) )
                points( sigma2S, gammaP, cex=scalefac*rib, pch=21, col='grey', bg=adjustcolor('blue',.3) )
                points( sigma2S, gammaP, cex=scalefac*sub.pelvic, pch=21, col='grey', bg=adjustcolor('grey',.3) )
            } )
}

