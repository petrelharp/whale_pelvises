#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(Matrix)

runfiles <- list.files(path="shape-likelihood-surface",pattern=".*RData",full.names=TRUE)
runs <- lapply( runfiles, function (x) { tmpenv <- environment(); load(x,envir=tmpenv); as.list(tmpenv) } )

pargrids <- lapply( runs, '[[', 'pargrid' )
grid1 <- ( sapply(pargrids,nrow) == 6400 )
pargrid1 <- pargrids[[ min(which(grid1)) ]]
pargrid2 <- pargrids[[ min(which(!grid1)) ]]
tmp1 <- cbind( pargrid1, kk=1:nrow(pargrid1) )
tmp2 <- cbind( pargrid2, kk=1:nrow(pargrid2) )
tmp <- merge(tmp1,tmp2,by=names(pargrid1),all=TRUE)
tmp <- tmp[do.call(order,tmp[1:3]),]
pargrid <- tmp[,1:3]
sigma2S.vals <- sort( unique( pargrid$sigma2S ) )
gammaP.vals <- sort( unique( pargrid$gammaP ) )

pelvic.sgrids1 <- do.call( cbind, do.call( c, lapply( runs[grid1], "[[", "pelvic.sgrids" ) ) )
pelvic.sgrids2 <- do.call( cbind, do.call( c, lapply( runs[!grid1], "[[", "pelvic.sgrids" ) ) )
pelvic.sgrids <- cbind( pelvic.sgrids1[tmp$kk.x,], pelvic.sgrids2[tmp$kk.y] )
sub.pelvic.sgrids1 <- do.call( cbind, do.call( c, lapply( runs[grid1], "[[", "sub.pelvic.sgrids" ) ) )
sub.pelvic.sgrids2 <- do.call( cbind, do.call( c, lapply( runs[!grid1], "[[", "sub.pelvic.sgrids" ) ) )
sub.pelvic.sgrids <- cbind( sub.pelvic.sgrids1[tmp$kk.x,], sub.pelvic.sgrids2[tmp$kk.y] )
rib.sgrids1 <- do.call( cbind, do.call( c, lapply( runs[grid1], "[[", "rib.sgrids" ) ) )
rib.sgrids2 <- do.call( cbind, do.call( c, lapply( runs[!grid1], "[[", "rib.sgrids" ) ) )
rib.sgrids <- cbind( rib.sgrids1[tmp$kk.x,], rib.sgrids2[tmp$kk.y] )


## normalize, within each
pelvic.normed <- exp(pelvic.sgrids-max(pelvic.sgrids,na.rm=TRUE)) 
pelvic.normed <- sweep( pelvic.normed, 2, colSums(pelvic.normed,na.rm=TRUE), "/" )
rib.normed <- exp(rib.sgrids-max(rib.sgrids,na.rm=TRUE)) 
rib.normed <- sweep( rib.normed, 2, colSums(rib.normed,na.rm=TRUE), "/" )
sub.pelvic.normed <- exp(sub.pelvic.sgrids-max(sub.pelvic.sgrids,na.rm=TRUE)) 
sub.pelvic.normed <- sweep( sub.pelvic.normed, 2, colSums(sub.pelvic.normed,na.rm=TRUE), "/" )

pargrid$pelvic <- rowMeans( pelvic.normed,na.rm=TRUE )
pargrid$rib <- rowMeans( rib.normed,na.rm=TRUE )
pargrid$sub.pelvic <- rowMeans( sub.pelvic.normed,na.rm=TRUE )
scalefac <- 100

##
# marginal likelihoods of ks
sapply( c("pelvic","sub.pelvic","rib"), function (x) tapply( pargrid[[x]]/sum(pargrid[[x]]), pargrid$ks, sum, na.rm=TRUE ) )

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
            tapply( pargrid[[type]], pargrid$gammaP, sum, na.rm=TRUE ) } )
plot( 0, type='n', xlab='gammaP', ylab='probability', xlim=range(pargrid$gammaP), ylim=range(gP.dists) )
fillcols <- adjustcolor(c('red','blue','grey'),.3)
for (k in 1:3) {
    polygon( c( min(gammaP.vals), gammaP.vals, max(gammaP.vals) ), c(0,gP.dists[[k]],0), col=fillcols[k] )
}
sS.dists <- lapply( c("pelvic","sub.pelvic","rib"), function (type) {
            tapply( pargrid[[type]], pargrid$sigma2S, sum, na.rm=TRUE ) } )
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
