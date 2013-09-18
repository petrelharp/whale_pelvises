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
        plot( shape_difference ~ treedist, data=shapediff, subset=(bone1==type), col=adjustcolor(c('black','red'),.2)[ifelse(sex1==sex2,1,2)], pch=20, main=type )
        ltmp <- loess( shape_difference ~ treedist, data=shapediff, subset=(bone1==type & sex1==sex2), span=.3 )
        lines( seq( 0, max(shapediff$treedist), length.out=100 ), predict( ltmp, newdata=data.frame(treedist=seq( 0, max(shapediff$treedist), length.out=100 )) ) )
        ltmp <- loess( shape_difference ~ treedist, data=shapediff, subset=(bone1==type & sex1!=sex2), span=.3 )
        lines( seq( 0, max(shapediff$treedist), length.out=100 ), predict( ltmp, newdata=data.frame(treedist=seq( 0, max(shapediff$treedist), length.out=100 )) ), col='red' )
    }
}

## sample-by-sample
# matrix of pelvis shape differences: exclude within-individual comparisons
pelvs <- with(shapediff, ( bone1 == "pelvic" ) & ( bone2 == "pelvic" ) & ( specimen1 != specimen2 ) )
pelvicdiff <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
pelvicdiff[ cbind(as.numeric(shapediff$specimen1)[pelvs],as.numeric(shapediff$specimen2)[pelvs]) ] <- shapediff$shape_difference[pelvs]
# matrix of rib shape differences: exclude within-individual comparisons
ribs <- with( shapediff, ( bone1 == "rib" ) & ( bone2 == "rib" ) & ( specimen1 != specimen2 ) )
ribdiff <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
ribdiff[ cbind(as.numeric(shapediff$specimen1)[ribs],as.numeric(shapediff$specimen2)[ribs]) ] <- shapediff$shape_difference[ribs]

###
# species-by-species
# rib
rib.speciesdiff <- with( subset(shapediff, bone1==bone2 & bone1=="rib" & ( specimen1 != specimen2 ) ), 
        tapply( shape_difference, list(species1,species2), sum, na.rm=TRUE ) )
rib.speciesdiff[is.na(rib.speciesdiff)] <- 0
rib.speciesdiff <- rib.speciesdiff + t(rib.speciesdiff) - diag(diag(rib.speciesdiff))
rib.denom <- with( subset(shapediff, bone1==bone2 & bone1=="rib" & ( specimen1 != specimen2 ) ), 
        tapply( shape_difference, list(species1,species2), length ) )
rib.denom[is.na(rib.denom)] <- 0
rib.denom <- rib.denom + t(rib.denom) - diag(diag(rib.denom))
rib.speciesdiff <- rib.speciesdiff / rib.denom
# pelvis
pelvic.speciesdiff <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic" & ( specimen1 != specimen2 ) ), 
        tapply( shape_difference, list(species1,species2), sum, na.rm=TRUE ) )
pelvic.speciesdiff[is.na(pelvic.speciesdiff)] <- 0
pelvic.speciesdiff <- pelvic.speciesdiff + t(pelvic.speciesdiff) - diag(diag(pelvic.speciesdiff))
pelvic.denom <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic" & ( specimen1 != specimen2 ) ), 
        tapply( shape_difference, list(species1,species2), length ) )
pelvic.denom[is.na(pelvic.denom)] <- 0
pelvic.denom <- pelvic.denom + t(pelvic.denom) - diag(diag(pelvic.denom))
pelvic.speciesdiff <- pelvic.speciesdiff / pelvic.denom


######
# get initial parameters
tdistmat <- with( shapediff, tapply( treedist, list(species1,species2), mean ) )
ut <- upper.tri(tdistmat,diag=TRUE)

nls.fits <- lapply( c(rib='rib',pelvic='pelvic'), function (type) {
        sdiff <- get(paste(type,".speciesdiff",sep=''))[ut]
        tdiff <- tdistmat[ut]
        nls( sdiff^2 ~ 2 * xi2P + ksig * tdiff, start=list(xi2P=.8,ksig=2) )
    } )

if (interactive()) {
    tvals <- seq(0,max(shapediff$treedist),length.out=100)
    layout(matrix(1:4,nrow=2))
    for (type in c('rib','pelvic')) {
        with(subset(shapediff,bone1==type), plot(treedist, shape_difference, cex=.25, pch=20, log='y' ) )
        with(subset(shapediff,bone1==type), lines( tvals, predict(loess(shape_difference ~ treedist, subset=(treedist>0)),newdata=data.frame(treedist=tvals)), col='red' ) )
        points( tdistmat[ut], get(paste(type,".speciesdiff",sep=''))[ut], col='red' )
        with( as.list(coef(nls.fits[[type]])), lines( tvals, sqrt( 2*xi2P + tvals*ksig ), col='green', lty=2 ) )
        with(subset(shapediff,bone1==type), plot(treedist, shape_difference, cex=.25, pch=20, xlim=c(0,.08), log='y' ) )
        with(subset(shapediff,bone1==type), lines( tvals, predict(loess(shape_difference ~ treedist, subset=(treedist>0)),newdata=data.frame(treedist=tvals)), col='red' ) )
        points( tdistmat[ut], get(paste(type,".speciesdiff",sep=''))[ut], col='red' )
        with( as.list(coef(nls.fits[[type]])), lines( tvals, sqrt( 2*xi2P + tvals*ksig ), col='green', lty=2 ) )
    }
}

initpar <- lapply( nls.fits, coef )

####
# setup for likelihood on the species tree
ut <- upper.tri(pelvic.speciesdiff)
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
        jx <- row(pelvic.speciesdiff)[ut][j]
        jy <- col(pelvic.speciesdiff)[ut][j]
        kx <- row(pelvic.speciesdiff)[ut][k]
        ky <- col(pelvic.speciesdiff)[ut][k]
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


#####
## just save everything
save( list=ls(), file="shape-brlens-stuff.RData" )

