#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(colorspace)
require(Matrix)

####
## read in data
load("tree-stuff.RData")
load("thedata-and-covmatrices.Rdata")
load("mcmc-setup.RData")
havedata <- !is.na(thedata)
datavec <- crossprod( projmatrix[havedata,], thedata[havedata] )   # true data
whales <- read.csv("whales.csv",header=TRUE)

# shapes
shapediff <- read.table("61_pairwise_shapes.out", header=TRUE)
### TEMPORARY: exclude problem sample
shapediff <- subset( shapediff, (specimen1 != "USNM_504107") & (specimen2 != "USNM_504107")  )
####
shapediff$specimen1 <- factor( shapediff$specimen1, levels=rownames(thedata) )
shapediff$specimen2 <- factor( shapediff$specimen2, levels=rownames(thedata) )
shapediff$species1 <- whales$species[match(shapediff$specimen1,whales$specimen)]
shapediff$species2 <- whales$species[match(shapediff$specimen2,whales$specimen)]
# restrict to same-bone comparisons and ones observed in the data
shapediff <- subset( shapediff, ( bone1==bone2 ) & ( specimen1 %in% rownames(thedata) ) & ( specimen2 %in% rownames(thedata) ) )

# want SQUARED distance
shapediff$sqdist <- shapediff$shape_difference^2

# remove females
no.females <- TRUE
if (no.females) {
    males <- whales$specimen[ whales$sex == "M" ]
    shapediff <- subset( shapediff, specimen1 %in% males & specimen2 %in% males )
}

# treedist
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
all.treedist <- treedist( adjtree )
rownames(all.treedist) <- colnames(all.treedist) <- rownames(thedata)
species.order <- match( levels(shapediff$species1), species.tree$tip.label  )
species.treedist <- treedist(species.tree)[ species.order, species.order ]
rownames(species.treedist) <- colnames(species.treedist) <- species.tree$tip.label[species.order]

# add variables
shapediff$treedist <- species.treedist[ as.matrix( shapediff[c("species1","species2")] ) ]
shapediff$sex1 <- whales$sex[ match(shapediff$specimen1,whales$specimen) ]
shapediff$sex2 <- whales$sex[ match(shapediff$specimen2,whales$specimen) ]
if (interactive() & !no.females) {
    layout(1:2)
    for (type in c('rib','pelvic')) {
        plot( shape_difference ~ jitter(treedist), data=shapediff, subset=(bone1==type), col=adjustcolor(c('blue','green','black'),.2)[ifelse(sex1==sex2,ifelse(sex1=='M',1,2),3)], pch=20, main=type, cex=.5 )
        ltmp <- loess( shape_difference ~ treedist, data=shapediff, subset=(bone1==type & sex1==sex2 & sex1=='M'), span=.3 )
        lines( seq( 0, max(shapediff$treedist), length.out=100 ), predict( ltmp, newdata=data.frame(treedist=seq( 0, max(shapediff$treedist), length.out=100 )) ), col='blue' )
        ltmp <- loess( shape_difference ~ treedist, data=shapediff, subset=(bone1==type & sex1==sex2 & sex1=='F'), span=.3 )
        lines( seq( 0, max(shapediff$treedist), length.out=100 ), predict( ltmp, newdata=data.frame(treedist=seq( 0, max(shapediff$treedist), length.out=100 )) ), col='green' )
        ltmp <- loess( shape_difference ~ treedist, data=shapediff, subset=(bone1==type & sex1!=sex2), span=.3 )
        lines( seq( 0, max(shapediff$treedist), length.out=100 ), predict( ltmp, newdata=data.frame(treedist=seq( 0, max(shapediff$treedist), length.out=100 )) ), col='black' )
    }
}

## sample-by-sample
# matrix of pelvis shape differences: exclude within-individual comparisons
pelvs <- with(shapediff, ( bone1 == "pelvic" ) & ( bone2 == "pelvic" ) & ( specimen1 != specimen2 ) )
pelvicdiffsq <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
pelvicdiffsq[ cbind(as.numeric(shapediff$specimen1)[pelvs],as.numeric(shapediff$specimen2)[pelvs]) ] <- shapediff$sqdist[pelvs]
# matrix of rib shape differences: exclude within-individual comparisons
ribs <- with( shapediff, ( bone1 == "rib" ) & ( bone2 == "rib" ) & ( specimen1 != specimen2 ) )
ribdiffsq <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
ribdiffsq[ cbind(as.numeric(shapediff$specimen1)[ribs],as.numeric(shapediff$specimen2)[ribs]) ] <- shapediff$sqdist[ribs]

###
# species-by-species
# rib
rib.speciesdiffsq <- with( subset(shapediff, bone1==bone2 & bone1=="rib" & ( specimen1 != specimen2 ) ), 
        tapply( sqdist, list(species1,species2), sum, na.rm=TRUE ) )
rib.speciesdiffsq[is.na(rib.speciesdiffsq)] <- 0
rib.speciesdiffsq <- rib.speciesdiffsq + t(rib.speciesdiffsq) - diag(diag(rib.speciesdiffsq))
rib.denom <- with( subset(shapediff, bone1==bone2 & bone1=="rib" & ( specimen1 != specimen2 ) ), 
        tapply( sqdist, list(species1,species2), length ) )
rib.denom[is.na(rib.denom)] <- 0
rib.denom <- rib.denom + t(rib.denom) - diag(diag(rib.denom))
rib.speciesdiffsq <- rib.speciesdiffsq / rib.denom
# pelvis
pelvic.speciesdiffsq <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic" & ( specimen1 != specimen2 ) ), 
        tapply( sqdist, list(species1,species2), sum, na.rm=TRUE ) )
pelvic.speciesdiffsq[is.na(pelvic.speciesdiffsq)] <- 0
pelvic.speciesdiffsq <- pelvic.speciesdiffsq + t(pelvic.speciesdiffsq) - diag(diag(pelvic.speciesdiffsq))
pelvic.denom <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic" & ( specimen1 != specimen2 ) ), 
        tapply( sqdist, list(species1,species2), length ) )
pelvic.denom[is.na(pelvic.denom)] <- 0
pelvic.denom <- pelvic.denom + t(pelvic.denom) - diag(diag(pelvic.denom))
pelvic.speciesdiffsq <- pelvic.speciesdiffsq / pelvic.denom


######
# get initial parameters
tdistmat <- with( shapediff, tapply( treedist, list(species1,species2), mean ) )
ut <- upper.tri(tdistmat,diag=TRUE)

nls.fits <- lapply( c(rib='rib',pelvic='pelvic'), function (type) {
        sdiff <- get(paste(type,".speciesdiffsq",sep=''))[ut]
        tdiff <- tdistmat[ut]
        nls( sdiff ~ 2 * kxi2P + ksig * tdiff, start=list(kxi2P=.8,ksig=2) )
    } )

optim.fits <- lapply( c(rib='rib',pelvic='pelvic'), function (type) {
        datavec <- subset( shapediff, bone1==bone2 & bone1==type )$sqdist
        treevec <- subset( shapediff, bone1==bone2 & bone1==type )$treedist
        llikfun <- function (par) {
            #  par = (xi2P, ks, sigma2P) 
            xi2P <- par[1]; ks <- par[2]; sigma2P <- par[3]
            brlens <- 2 * xi2P + sigma2P * treevec
            retval <- (-1) * sum( dchisq( datavec/brlens, df=ks, log=TRUE ) ) + sum( log(brlens) )
            if (!is.finite(retval)) { browser() } else { return(retval) }
        }
        ans <- optim( par=c( xi2P=coef(nls.fits[[type]])['kxi2P']/20, ks=20, sigma2P=coef(nls.fits[[type]])['ksig']/20 ), fn=llikfun, lower=c(0.01,0.1,0.1), method="L-BFGS-B", control=list(parscale=c(.1,.1,.1),fnscale=1e4) )
        names(ans$par) <- c('xi2P','ks','sigma2P')
        ans
    } )

if (interactive()) {
    tvals <- seq(0,max(shapediff$treedist),length.out=100)
    layout(matrix(1:4,nrow=2))
    for (type in c('rib','pelvic')) {
        with(subset(shapediff,bone1==type), plot(treedist, sqdist, cex=.25, pch=20, log='y', main=type, ylab='squared distance', xlab='tree distance' ) )
        with(subset(shapediff,bone1==type), lines( tvals, predict(loess(sqdist ~ treedist, subset=(treedist>0)),newdata=data.frame(treedist=tvals)), col='red' ) )
        points( tdistmat[ut], get(paste(type,".speciesdiffsq",sep=''))[ut], col='red' )
        with( as.list(coef(nls.fits[[type]])), lines( tvals, ( 2*kxi2P + tvals*ksig ), col='green', lty=2 ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, ( 2*xi2P + tvals ) * ks * sigma2P, col='purple' ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, qchisq( p=.25, df=ks ) * ( 2*xi2P + tvals ), col='purple', lty=3 ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, qchisq( p=.75, df=ks ) * ( 2*xi2P + tvals ), col='purple', lty=3 ) )
        with(subset(shapediff,bone1==type), plot(treedist, sqdist, cex=.25, pch=20, xlim=c(0,.08), log='y', main=type, ylab='squared distance', xlab='tree distance' ) )
        with(subset(shapediff,bone1==type), lines( tvals, predict(loess(sqdist ~ treedist, subset=(treedist>0)),newdata=data.frame(treedist=tvals)), col='red' ) )
        points( tdistmat[ut], get(paste(type,".speciesdiffsq",sep=''))[ut], col='red' )
        with( as.list(coef(nls.fits[[type]])), lines( tvals, ( 2*kxi2P + tvals*ksig ), col='green', lty=2 ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, ( 2*xi2P + tvals ) * ks * sigma2P, col='purple' ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, qchisq( p=.25, df=ks ) * ( 2*xi2P + tvals ), col='purple', lty=3 ) )
        with( as.list( optim.fits[[type]]$par ), lines( tvals, qchisq( p=.75, df=ks ) * ( 2*xi2P + tvals ), col='purple', lty=3 ) )
    }
}

initpar <- lapply( optim.fits, "[[", 'par' )

# $rib:  xi2P         ks    sigma2P 
# 0.05750125 5.93839502 0.43632759 
# $pelvic:  xi2P         ks    sigma2P 
# 0.09504044 8.88862550 0.97149192 #
# check xi2P is reasonable
(1/2) * with( subset(shapediff,species1==species2), tapply( sqdist, list(bone1,ifelse(specimen1==specimen2,'same','diff')), mean ) )
#  ... pretty close... within-individual variance is about 1/3.5 times this.

####
# setup for likelihood on the species tree
#   shared.paths is matrix of shared path lengths for pairs of comparisons
# also, tip length parameters:
#   *.shared.indivs gives the sum of shared individuals across all comparisons
#   *.shared.samples gives the sum of shared samples across all comparisons
#   *.shared.samples.indivs gives the sum of the product of shared samples and shared indivs, across all comparisons
#   *.shared.indivs.sq gives the sum of shared individuals across all comparisons, squared
#   *.shared.samples.sq gives the sum of shared samples across all comparisons, squared
# THEN
#   covariances between the means of all species-species comparisons is,
#     with y = shared.indivs * xi2i
#          z = shared.samples * xi2s
#     and y.z = shared.samples.indivs * xi2i * xi2s
#     and y.sq = shared.indivs.sq * xi2i^2
#          z = shared.samples.sq * xi2s^2
#     shared.paths^2 + 2 * shared.paths * (y+z) + 2 * y.z + y.sq + z.sq
ut <- upper.tri(pelvic.speciesdiffsq)
## shared.paths
# construct covariance matrix
sptree <- drop.tip( species.tree, setdiff( species.tree$tip.label, levels(shapediff$species1) ) )
spdesc <- get.descendants(sptree)
sp.edge.indices <- sptree$edge[,2] # associate each edge with the downstream node
sp.tip.edges <- match( 1:Ntip(sptree), sp.edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
sp.rootnode <- which.max( rowSums(spdesc) )  # number of tips each edge contributes to
shared.paths <- matrix( 0, nrow=choose(Ntip(sptree),2), ncol=choose(Ntip(sptree),2) ) # in order as upper.tri
sput <- upper.tri(shared.paths,diag=TRUE)
# matrix so that elements of upper.tri(shared.paths) are sp.mapping %*% c(edge.length,0)
if (!file.exists("spmapping.RData")) {  # takes a minute
    sp.mapping <- matrix( 0, nrow=sum(sput), ncol=Nedge(sptree) )
    m <- 1
    for (j in 1:nrow(shared.paths))  for (k in 1:j) {
        jx <- row(pelvic.speciesdiffsq)[ut][j]
        jy <- col(pelvic.speciesdiffsq)[ut][j]
        kx <- row(pelvic.speciesdiffsq)[ut][k]
        ky <- col(pelvic.speciesdiffsq)[ut][k]
        jxdesc <- spdesc[,jx] # which nodes are along the path from the root to each of these tips
        jydesc <- spdesc[,jy] 
        kxdesc <- spdesc[,kx] 
        kydesc <- spdesc[,ky] 
        #  which nodes are on both of the paths between the two pairs of tips
        contributes.nodes <- which( ( ( jxdesc | jydesc ) & ! ( jxdesc & jydesc ) ) & ( ( kxdesc | kydesc ) & ! ( kxdesc & kydesc ) ) )
        contributes.edges <- sp.edge.indices %in% contributes.nodes  # translate to edges
        sp.mapping[m,] <- contributes.edges
        m <- m+1
    }
    stopifnot( m == sum(sput) + 1 )
    spmap.nonz <- ( rowSums(sp.mapping) > 0 )
    sp.mapping <- Matrix( sp.mapping[spmap.nonz,] )
    save(sp.mapping, spmap.nonz, file="spmapping.RData")
} else {
    load("spmapping.RData")
}
shared.paths[ upper.tri(shared.paths,diag=TRUE) ][spmap.nonz] <- as.vector( sp.mapping %*% sptree$edge.length )  ### update like this, but with @x
shared.paths[lower.tri(shared.paths)] <- t(shared.paths)[lower.tri(shared.paths)]
shared.paths <- Matrix( shared.paths )
stopifnot( class(shared.paths) == "dsyMatrix"  & shared.paths@uplo == "U" )  # if so, changing upper tri also changes lower tri
###
# *.shared.{indivs,samples}.*
rownames.ut <- rownames(pelvic.speciesdiffsq)[row(pelvic.speciesdiffsq)[ut]]
colnames.ut <- colnames(pelvic.speciesdiffsq)[col(pelvic.speciesdiffsq)[ut]]
both.specimens <- intersect( with(subset(shapediff,bone1=='rib'),c(specimen1,specimen2)), with(subset(shapediff,bone1=='pelvic'),c(specimen1,specimen2)) )
if (!file.exists("shared-num-bones.RData")) {
    for (type in c('rib','pelvic')) {
        sharenames <- c( "shared.indivs", "shared.samples", "shared.samples.indivs", "shared.samples.sq", "shared.indivs.sq" )
        for (x in sharenames) { assign(x,matrix( 0, nrow=choose(Ntip(sptree),2), ncol=choose(Ntip(sptree),2) ) ) }
        m <- 1
        for (j in 1:nrow(shared.paths))  for (k in 1:j) {
            jx <- rownames.ut[j]
            jy <- colnames.ut[j]
            kx <- rownames.ut[k]
            ky <- colnames.ut[k]
            if ( length( intersect( c(jx,jy), c(kx,ky) ) ) > 0 ) {
                j.obs <- subset( shapediff, ( ( species1==jx & species2==jy ) | ( species2==jx & species1==jy ) ) & ( bone1 == type ) )
                k.obs <- subset( shapediff, ( ( species1==kx & species2==ky ) | ( species2==kx & species1==ky ) ) & ( bone1 == type ) )
                jk.denom <- ( nrow(j.obs) * nrow(k.obs) ) / ( 1 + (j==k) )
                if ( jk.denom > 0 ) {
                    shindivs <- outer( 1:nrow(j.obs), 1:nrow(k.obs), function (jj,kk) { 
                                ( j.obs$specimen1[jj] == k.obs$specimen1[kk] ) + 
                                ( j.obs$specimen1[jj] == k.obs$specimen2[kk] ) + 
                                ( j.obs$specimen2[jj] == k.obs$specimen1[kk] ) + 
                                ( j.obs$specimen2[jj] == k.obs$specimen2[kk] ) } )
                    shsamples <- outer( 1:nrow(j.obs), 1:nrow(k.obs), function (jj,kk) { 
                                ( j.obs$specimen1[jj] == k.obs$specimen1[kk] ) & ( j.obs$side1[jj] == k.obs$side1[kk] ) + 
                                ( j.obs$specimen1[jj] == k.obs$specimen2[kk] ) & ( j.obs$side1[jj] == k.obs$side2[kk] ) + 
                                ( j.obs$specimen2[jj] == k.obs$specimen1[kk] ) & ( j.obs$side2[jj] == k.obs$side1[kk] ) + 
                                ( j.obs$specimen2[jj] == k.obs$specimen2[kk] ) & ( j.obs$side2[jj] == k.obs$side2[kk] ) } )
                    shared.indivs[j,k] <- sum( shindivs ) / jk.denom^2
                    shared.samples[j,k] <- sum( shsamples ) / jk.denom^2
                    shared.samples.indivs[j,k] <- sum( shindivs*shsamples ) / jk.denom^2
                    shared.indivs.sq[j,k] <- sum( shindivs^2 ) / jk.denom^2
                    shared.samples.sq[j,k] <- sum( shsamples^2 ) / jk.denom^2
                }
            }
            m <- m+1
        }
        for (x in sharenames) { 
            z <- get(x)
            assign( paste(type,x,sep='.'), Matrix(z + t(z) - diag(diag(z))) )
        }
    }
    save( list=outer( c('rib','pelvic'), sharenames, paste, sep='.' ), file="shared-num-bones.RData" )
} else {
    load("shared-num-bones.RData")
}


#####
## just save everything
save( list=ls(), file="shape-brlens-stuff.RData" )

