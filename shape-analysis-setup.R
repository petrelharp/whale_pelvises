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

# parameter estimates from here
load("results.RData")
estpar <- unlist( pars["estpar",] )

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

# treedist
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
all.treedist <- treedist( adjtree )
rownames(all.treedist) <- colnames(all.treedist) <- rownames(thedata)
species.order <- match( levels(shapediff$species1), species.tree$tip.label  )
species.treedist <- treedist(species.tree)[ species.order, species.order ]
rownames(species.treedist) <- colnames(species.treedist) <- species.tree$tip.label[species.order]


####
# get length-adjusted testes size on the branches:
#  (posterior mean of testes size) - (predicted testes size given posterior mean of length)
fullmat <- make.fullmat( estpar )
# posterior mean trait values relative to the root at all nodes:
get.postmean <- function (usethese) {
    themean <- sweep( thedata, 2, unlist(phylomeans), "-" )
    themean[!usethese] <- fullmat[!usethese,usethese] %*% solve( fullmat[usethese,usethese], themean[usethese] )
    dimnames(themean) <- dimnames(thedata)
    themean <- sweep( themean, 2, unlist(phylomeans), "+" )
    return( themean )
}
fullmean <- get.postmean(havedata)
islength <- ( col(thedata) == which(colnames(thedata)=="bodylength") & havedata )
istestes <- ( col(thedata) == which(colnames(thedata)=="actual_testes_mass_max") )
lengthmean <- get.postmean( islength )
# testes:
fullmean <- cbind( fullmean, scaled.testes=(fullmean[,"actual_testes_mass_max"] - lengthmean[,"actual_testes_mass_max"]) )
# plot( lengthmean[,"actual_testes_mass_max"], fullmean[,'scaled.testes'] ) ## looks good
# associate mean value of adjacent nodes to each edge
edge.values <- t( apply( tree$edge, 1, function (kk) { colMeans(fullmean[kk,]) } ) )
colnames(edge.values) <- colnames(fullmean)
edge.testes <- ( edge.values[,"scaled.testes"] )

sample.values <- function (dothese=(istestes & ! havedata)) {
    withthese <- ( ! dothese ) & havedata
    spar <- samples[ sample(1:nrow(samples),1), ]
    lengthmean <- get.postmean( islength )
    sampled.cov <- fullmat[dothese,dothese] - fullmat[dothese,withthese] %*% solve( fullmat[withthese,withthese], fullmat[withthese,dothese] ) 
    schol <- chol(sampled.cov,pivot=TRUE) 
    sampled.values <- get.postmean(withthese) 
    sampled.values[dothese] <- sampled.values[dothese] + t( schol[,order(attr(schol,"pivot"))] ) %*% rnorm(sum(dothese))
    sampled.values <- cbind( sampled.values, scaled.testes=(sampled.values[,"actual_testes_mass_max"] - lengthmean[,"actual_testes_mass_max"]) )
    edge.values <- t( apply( tree$edge, 1, function (kk) { colMeans(sampled.values[kk,]) } ) )
    colnames(edge.values) <- colnames(sampled.values)
    return( edge.values[,"scaled.testes"] )
}

sampled.edge.testes <- t( replicate(1000, sample.values() ) )

# matplot(t(sweep(sampled.edge.testes,2,edge.testes,"-")),type='l')


# and testes-weighted relative time in the tree:
internal.lengths <- tree$edge.length
internal.lengths[tip.edges] <- (.05/3.16)^2
testes.treedist <- treedist( tree, edge.length=scale(edge.testes * internal.lengths) )
rownames(testes.treedist) <- colnames(testes.treedist) <- rownames(thedata)
testes.spdist <- testes.treedist[match(levels(shapediff$species1),rownames(thedata)),match(levels(shapediff$species1),rownames(thedata))]

# visualization

if (interactive()) {
    layout(t(1:2))
    plot(adjtree,edge.color=diverge_hcl(64)[cut(edge.testes,breaks=seq((-1.05)*max(abs(edge.testes)),1.05*max(abs(edge.testes)),length.out=65))],main='relative testes size',edge.width=3,tip.color=ifelse(is.na(thedata[,"left.rib"])&is.na(thedata[,"right.rib"]),'red','black')[1:Ntip(adjtree)],cex=.5)
    plot(adjtree,edge.color=diverge_hsv(64)[cut(fullmean[edge.indices,'bodylength'],64)],main='body length',edge.width=3,cex=.5,tip.color=ifelse(is.na(thedata[,"left.rib"])&is.na(thedata[,"right.rib"]),'red','black')[1:Ntip(adjtree)])
}


######
# BACK-OF-THE-ENVELOPE:

###
# all the data
pelvics <- with(shapediff, (bone1 == bone2) & (bone1 == "pelvic") & (specimen1 != specimen2) )
pelvic.specimens <- with( subset(shapediff,pelvics), unique( c(as.character(specimen1),as.character(specimen2)) ) )
ribs <- with(shapediff, (bone1 == bone2) & (bone1 == "rib") & (specimen1 != specimen2) )
rib.specimens <- with( subset(shapediff,ribs), unique( c(as.character(specimen1),as.character(specimen2)) ) )
both.specimens <- intersect( pelvic.specimens, rib.specimens )
have.both <- with(shapediff, specimen1 %in% both.specimens &  specimen2 %in% both.specimens )
havecols <- ifelse( have.both, 'black', 'red' )

tdist <- all.treedist[ as.matrix( shapediff[,c("specimen1","specimen2")] ) ]
testdist <- testes.treedist[ as.matrix( shapediff[,c("specimen1","specimen2")] ) ]
pelvic.lm <- with( shapediff, lm( shape_difference ~ tdist, subset=pelvics ) )
rib.lm <- with( shapediff, lm( shape_difference ~ tdist, subset=ribs ) )
pelvic.bi.lm <- with( shapediff, lm( shape_difference ~ tdist + testdist, subset=pelvics ) )
rib.bi.lm <- with( shapediff, lm( shape_difference ~ tdist + testdist, subset=ribs ) )

if (interactive()) {
    layout(t(1:2))
    with(shapediff, plot( shape_difference ~ tdist, pch=20, cex=.5, col=adjustcolor(havecols,.2), subset=pelvics, main='pelvics' ) )
    points( fitted(pelvic.bi.lm) ~ tdist[pelvics], cex=.5, col=adjustcolor('green',.2) )
    abline( coef(pelvic.lm), col='red' )
    with(shapediff, plot( shape_difference ~ tdist, pch=20, cex=.5, col=adjustcolor(havecols,.2), subset=ribs, main='ribs' ) )
    points( fitted(rib.bi.lm) ~ tdist[ribs], cex=.5, col=adjustcolor('green',.2) )
    abline( coef(rib.lm), col='red' )
}


##
# species means

# predict shape difference from time in the tree:
ut <- upper.tri(pelvic.speciesdiff)
have.ribs <- ( !is.na( rib.speciesdiff[ut] ) )
have.pelvics <- ( !is.na( pelvic.speciesdiff[ut] ) )
have.both <- have.ribs & have.pelvics
pelvic.lm <- lm( pelvic.speciesdiff[ut] ~ species.treedist[ut] )
bi.pelvic.lm <- lm( pelvic.speciesdiff[ut] ~ species.treedist[ut] * testes.spdist[ut] )
rib.lm <- lm( rib.speciesdiff[ut] ~ species.treedist[ut] )
bi.rib.lm <- lm( rib.speciesdiff[ut] ~ species.treedist[ut] * testes.spdist[ut] )

# restrict to full observations
sub.pelvic.lm <- lm( pelvic.speciesdiff[ut] ~ species.treedist[ut] , subset=have.both )
sub.bi.pelvic.lm <- lm( pelvic.speciesdiff[ut] ~ species.treedist[ut] + testes.spdist[ut], subset=have.both )
sub.rib.lm <- lm( rib.speciesdiff[ut] ~ species.treedist[ut] , subset=have.both )
sub.bi.rib.lm <- lm( rib.speciesdiff[ut] ~ species.treedist[ut] + testes.spdist[ut], subset=have.both )

anova( bi.rib.lm, rib.lm )
anova( bi.pelvic.lm, pelvic.lm )

##########
# compare residuals to testes-weighted distance
rib.resids <- fitted(bi.rib.lm) - fitted(rib.lm)
pelvic.resids <- fitted(bi.pelvic.lm) - fitted(pelvic.lm)

resid.rib.lm <- lm( resid(rib.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.ribs] )
resid.pelvic.lm <- lm( resid(pelvic.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.pelvics] )
# restrict to full observations
sub.resid.rib.lm <- lm( resid(sub.rib.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.both] )
sub.resid.pelvic.lm <- lm( resid(sub.pelvic.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.both] )

if (interactive()) {
    # look at these:
    pdf( file="ribs-correlate.pdf", width=8, height=5, pointsize=10 )
    layout(matrix(1:4,nrow=2))
    cols <- ifelse(have.both,"black","red")
    ##
    # shape-tree-diffs
    plot( species.treedist[ut], rib.speciesdiff[ut], xlab="tree distance", ylab="shape difference", main="rib", col=cols[have.ribs], xlim=range(species.treedist) )
    abline(coef(rib.lm),col='red')
    abline(coef(sub.rib.lm))
    points( species.treedist[ut][have.ribs], fitted(bi.rib.lm), col='green', pch=20 )
    plot( species.treedist[ut], pelvic.speciesdiff[ut], xlab="tree distance", ylab="shape difference", main="pelvic", col=cols[have.pelvics], xlim=range(species.treedist) )
    # identify( species.treedist[ut], pelvic.speciesdiff[ut], labels=outer(rownames(pelvic.speciesdiff),colnames(pelvic.speciesdiff),paste)[ut] )
    abline(coef(pelvic.lm),col='red')
    abline(coef(sub.pelvic.lm))
    points( species.treedist[ut][have.pelvics], fitted(bi.pelvic.lm), col='green', pch=20 )
    # shape-testes-tree-diffs
    plot( testes.spdist[ut], rib.speciesdiff[ut], xlab="testes-weighted tree distance", ylab="shape difference", main="rib", col=cols[have.ribs], xlim=range(testes.spdist) )
    points( testes.spdist[ut][have.ribs], fitted(bi.rib.lm), col='green', pch=20 )
    plot( testes.spdist[ut], pelvic.speciesdiff[ut], xlab="testes-weighted tree distance", ylab="shape difference", main="pelvic", col=cols[have.pelvics], xlim=range(testes.spdist) )
    points( testes.spdist[ut][have.pelvics], fitted(bi.pelvic.lm), col='green', pch=20 )
    ##
    # resids-testes-distance
    plot( testes.spdist[ut][have.ribs], resid(rib.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="ribs", col=cols[have.ribs], xlim=range( testes.spdist ) )
    abline(coef(resid.rib.lm),col='red')
    plot( testes.spdist[ut][have.pelvics], resid(pelvic.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="pelvic, residuals of all data", col=cols[have.pelvics], xlim=range( testes.spdist ) )
    abline(coef(resid.pelvic.lm),col='red')
    # resids-testes-distance, all observations
    plot( testes.spdist[ut][have.both], resid(sub.rib.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="ribs", col=cols[have.both], xlim=range( testes.spdist ) )
    abline(coef(sub.resid.rib.lm))
    plot( testes.spdist[ut][have.both], resid(sub.pelvic.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="pelvic, residuals for complete obs", col=cols[have.both], xlim=range( testes.spdist ) )
    abline(coef(sub.resid.pelvic.lm))
    dev.off()
    ### Hm.  Looks good?
}



#####
# quick-and-dirty parameters
# within-species variance
tip.ribvar <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), tapply( shape_difference, species1, mean, na.rm=TRUE ) )
tip.pelvicvar <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), tapply( shape_difference, species1, mean, na.rm=TRUE ) )
# estimate of xi2P,R
est.xi2P <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), var( shape_difference, na.rm=TRUE ) / mean( shape_difference, na.rm=TRUE ) ) / 2
est.xi2R <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), var( shape_difference, na.rm=TRUE ) / mean( shape_difference, na.rm=TRUE ) ) / 2
# estimate of ksP,R
ksP <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), 2 * mean( shape_difference, na.rm=TRUE )^2 / var( shape_difference, na.rm=TRUE ) )
ksR <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), 2 * mean( shape_difference, na.rm=TRUE )^2 / var( shape_difference, na.rm=TRUE ) )

# check:
if (interactive()) {
    tmp <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), hist( shape_difference, main="pelvic", breaks=40 ) )
    lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2P, df=ksP )) )
    tmp <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), hist( shape_difference, main='rib', breaks=40 ) )
    lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2R, df=ksR )) )
    lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2R*(ksP/ksR), df=ksP )), col='red' )  # looks better
}

# estimates
est.ks <- ksP
est.xi2R <- est.xi2R * ksP/ksR

# sigma2S : slope of relation between shape diff and tree dist, divided by ks.
est.sigma2S <- (coef(pelvic.lm)/ksP)[2]
# gammaP : slope of relation to residual
est.gammaP <- (coef(resid.pelvic.lm)/ksP)[2]

initpar <- c( 
        ks = est.ks,
        sigma2S = est.sigma2S,
        gammaP = est.gammaP,
        xi2P = est.xi2P
    )
names(initpar) <- c('ks','sigma2S','gammaP','xi2P')


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

# get edge.testes on the species tree
## grr: map between the trees
descendants <- get.descendants(tree)
colnames( descendants ) <- rownames( descendants ) <- c( tree$tip.label, paste("NA",1:Nnode(tree),sep='.') )
tmp <- t( descendants[edge.indices,match(sptree$tip.label,tree$tip.label)] ) # edge-by-species-tip
spdesc <- get.descendants(sptree)
colnames( spdesc ) <- rownames( spdesc ) <- c( sptree$tip.label, paste("NA",1:Nnode(sptree),sep='.') )
sptmp <- t( spdesc[sp.edge.indices,1:Ntip(sptree)] )
## tree.translate[k] gives the index of the corresponding edge in tree
tree.translate <- lapply( 1:Nedge(sptree), function (k) {
        ematch <- which( apply( tmp == sptmp[,k], 2, all ) )
        return( ematch[ ! names(ematch) %in% sptree$tip.label ] )
    } )
stopifnot( all( sapply(tree.translate,length) == 1 ) )
tree.translate <- unlist(tree.translate)
## ok finally
sp.edge.testes <- edge.testes[ tree.translate ]
sampled.sp.edge.testes <- sampled.edge.testes[ , tree.translate ]
# check
#   layout(1:2)
#   plot(adjtree, edge.color=ifelse(edge.testes>0,'red','blue') )
#   plot(sptree, edge.color=ifelse(sp.edge.testes>0,'red','blue') )


#####
## just save everything
save( list=ls(), file="shape-stuff.RData" )


if (FALSE) {

########
# add tips?
rib.lengths <- tree$edge.length
rib.lengths[ tip.edges ] <- xi2R
pelvic.lengths <- tree$edge.length
pelvic.lengths[ tip.edges ] <- xi2P
rib.treemat <- treedist( tree, edge.length=rib.lengths )
pelvic.treemat <- treedist( tree, edge.length=pelvic.lengths )

layout(matrix(1:4,nrow=2))
plot( as.vector(rib.treemat), as.vector(ribdiff), ylab="rib differences" )
plot( as.vector(pelvic.treemat), as.vector(pelvicdiff), ylab="pelvis differences" )

}
