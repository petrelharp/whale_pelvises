#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)
require(colorspace)

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
species.order <- match( levels(shapediff$species1), species.tree$tip.label  )
species.treedist <- treedist( species.tree )[ species.order, species.order ]
rownames(species.treedist) <- colnames(species.treedist) <- species.tree$tip.label[species.order]


####
# get length-adjusted testes size on the branches:
#  (posterior mean of testes size) - (predicted testes size given posterior mean of length)
fullmat <- make.fullmat( estpar )
# posterior mean trait values relative to the root at all nodes:
fullmean <- sweep( thedata, 2, unlist(phylomeans), "-" )
fullmean[!havedata] <- fullmat[!havedata,havedata] %*% solve( fullmat[havedata,havedata], fullmean[havedata] )
dimnames(fullmean) <- dimnames(thedata)
fullmean <- sweep( fullmean, 2, unlist(phylomeans), "+" )
# predicted values given only length:
lengthmean <- sweep( thedata, 2, unlist(phylomeans), "-" )
islength <- ( col(thedata) == which(colnames(thedata)=="bodylength") & havedata )
lengthmean[!islength] <- fullmat[!islength,islength] %*% solve( fullmat[islength,islength], lengthmean[islength] )
dimnames(lengthmean) <- dimnames(thedata)
lengthmean <- sweep( lengthmean, 2, unlist(phylomeans), "+" )
# testes:
fullmean <- cbind( fullmean, scaled.testes=(fullmean[,"actual_testes_mass_max"] - lengthmean[,"actual_testes_mass_max"]) )
# plot( lengthmean[,"actual_testes_mass_max"], fullmean[,'scaled.testes'] ) ## looks good
# associate mean value of adjacent nodes to each edge
edge.values <- t( apply( tree$edge, 1, function (kk) { colMeans(fullmean[kk,]) } ) )
colnames(edge.values) <- colnames(fullmean)
edge.testes <- ( edge.values[,"scaled.testes"] * internal.lengths )

# and testes-weighted relative time in the tree:
testes.treedist <- treedist( tree, edge.length=scale(edge.testes) )
rownames(testes.treedist) <- colnames(testes.treedist) <- rownames(thedata)
testes.spdist <- testes.treedist[match(levels(shapediff$species1),rownames(thedata)),match(levels(shapediff$species1),rownames(thedata))]

# visualization
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
#
if (interactive()) {
    layout(t(1:2))
    plot(adjtree,edge.color=diverge_hcl(64)[cut(edge.testes,breaks=seq((-1.05)*max(abs(edge.testes)),1.05*max(abs(edge.testes)),length.out=65))],main='relative testes size',edge.width=3,tip.color=ifelse(is.na(thedata[,"left.rib"])&is.na(thedata[,"right.rib"]),'red','black')[1:Ntip(adjtree)],cex=.5)
    plot(adjtree,edge.color=diverge_hsv(64)[cut(fullmean[edge.indices,'bodylength'],64)],main='body length',edge.width=3,cex=.5,tip.color=ifelse(is.na(thedata[,"left.rib"])&is.na(thedata[,"right.rib"]),'red','black')[1:Ntip(adjtree)])
}


######
# BACK-OF-THE-ENVELOPE:

# predict shape difference from time in the tree:
have.ribs <- ( !is.na( rib.speciesdiff[upper.tri(rib.speciesdiff)] ) )
rib.lm <- ( lm( rib.speciesdiff[upper.tri(rib.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)], subset=have.ribs ) )
have.pelvics <- ( !is.na( pelvic.speciesdiff[upper.tri(pelvic.speciesdiff)] ) )
pelvic.lm <- ( lm( pelvic.speciesdiff[upper.tri(pelvic.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)], subset=have.pelvics ) )
# restrict to full observations
have.both <- ( !is.na( rib.speciesdiff[upper.tri(rib.speciesdiff)] ) ) & ( !is.na( pelvic.speciesdiff[upper.tri(pelvic.speciesdiff)] ) ) 
sub.rib.lm <- ( lm( rib.speciesdiff[upper.tri(rib.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)], subset=have.both ) )
sub.pelvic.lm <- ( lm( pelvic.speciesdiff[upper.tri(pelvic.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)], subset=have.both ) )

##########
# compare residuals to testes-weighted distance
resid.rib.lm <- lm( resid(rib.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.ribs] )
resid.pelvic.lm <- lm( resid(pelvic.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.pelvics] )
# restrict to full observations
sub.resid.rib.lm <- lm( resid(sub.rib.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.both] )
sub.resid.pelvic.lm <- lm( resid(sub.pelvic.lm) ~  testes.spdist[upper.tri(testes.spdist)][have.both] )

# look at these:
pdf( file="ribs-correlate.pdf", width=8, height=5, pointsize=10 )
layout(matrix(1:4,nrow=2))
cols <- ifelse(have.both,"black","red")
##
# shape-tree-diffs
plot( species.treedist[upper.tri(species.treedist)], rib.speciesdiff[upper.tri(species.treedist)], xlab="tree distance", ylab="shape difference", main="rib", col=cols[have.ribs], xlim=range(species.treedist) )
abline(coef(rib.lm),col='red')
abline(coef(sub.rib.lm))
plot( species.treedist[upper.tri(species.treedist)], pelvic.speciesdiff[upper.tri(species.treedist)], xlab="tree distance", ylab="shape difference", main="pelvic", col=cols[have.pelvics], xlim=range(species.treedist) )
identify( species.treedist[upper.tri(species.treedist)], pelvic.speciesdiff[upper.tri(species.treedist)], labels=outer(rownames(pelvic.speciesdiff),colnames(pelvic.speciesdiff),paste)[upper.tri(species.treedist)] )
abline(coef(pelvic.lm),col='red')
abline(coef(sub.pelvic.lm))
# shape-testes-tree-diffs
plot( testes.spdist[upper.tri(testes.spdist)], rib.speciesdiff[upper.tri(testes.spdist)], xlab="testes-weighted tree distance", ylab="shape difference", main="rib", col=cols[have.ribs], xlim=range(testes.spdist) )
plot( testes.spdist[upper.tri(testes.spdist)], pelvic.speciesdiff[upper.tri(testes.spdist)], xlab="testes-weighted tree distance", ylab="shape difference", main="pelvic", col=cols[have.pelvics], xlim=range(testes.spdist) )
##
# resids-testes-distance
plot( testes.spdist[upper.tri(testes.spdist)][have.ribs], resid(rib.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="ribs", col=cols[have.ribs], xlim=range( testes.spdist ) )
abline(coef(resid.rib.lm),col='red')
plot( testes.spdist[upper.tri(testes.spdist)][have.pelvics], resid(pelvic.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="pelvic, residuals of all data", col=cols[have.pelvics], xlim=range( testes.spdist ) )
abline(coef(resid.pelvic.lm),col='red')
# resids-testes-distance, all observations
plot( testes.spdist[upper.tri(testes.spdist)][have.both], resid(sub.rib.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="ribs", col=cols[have.both], xlim=range( testes.spdist ) )
abline(coef(sub.resid.rib.lm))
plot( testes.spdist[upper.tri(testes.spdist)][have.both], resid(sub.pelvic.lm), xlab="testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="pelvic, residuals for complete obs", col=cols[have.both], xlim=range( testes.spdist ) )
abline(coef(sub.resid.pelvic.lm))
dev.off()
### Hm.  Looks good?



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
tmp <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), hist( shape_difference, main="pelvic", breaks=40 ) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2P, df=ksP )) )
tmp <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), hist( shape_difference, main='rib', breaks=40 ) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2R, df=ksR )) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/est.xi2R*(ksP/ksR), df=ksP )), col='red' )  # looks better

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

#####
## just save everything
save( list=ls(), file="shape-stuff.RData" )
