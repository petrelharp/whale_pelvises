#/usr/bin/R --vanilla
source("correlated-traits-fns.R")
require(ape)

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
shapediff <- subset( shapediff, ( bone1==bone2 ) & ( specimen1 %in% rownames(thedata) ) & ( specimen2 %in% rownames(thedata) ) )
shapediff$species1 <- whales$species[match(shapediff$specimen1,whales$specimen)]
shapediff$species2 <- whales$species[match(shapediff$specimen2,whales$specimen)]
### TEMPORARY
shapediff <- subset( shapediff, (specimen1 != "USNM_504107") & (specimen2 != "USNM_504107") & ( side1 == "left") & ( side2 == "left" ) )
####
shapediff$specimen1 <- factor( shapediff$specimen1, levels=rownames(thedata) )
shapediff$specimen2 <- factor( shapediff$specimen2, levels=rownames(thedata) )
# matrix of pelvis shape differences
pelvs <- ( shapediff$bone1 == "pelvic" ) & ( shapediff$bone2 == "pelvic" )
pelvicdiff <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
pelvicdiff[ cbind(as.numeric(shapediff$specimen1)[pelvs],as.numeric(shapediff$specimen2)[pelvs]) ] <- shapediff$shape_difference[pelvs]
# matrix of rib shape differences
ribs <- ( shapediff$bone1 == "rib" ) & ( shapediff$bone2 == "rib" )
ribdiff <- matrix(NA, nrow=nrow(thedata), ncol=nrow(thedata) )
ribdiff[ cbind(as.numeric(shapediff$specimen1)[ribs],as.numeric(shapediff$specimen2)[ribs]) ] <- shapediff$shape_difference[ribs]

# by species
rib.speciesdiff <- with( subset(shapediff, bone1==bone2 & bone1=="rib"), tapply( shape_difference, list(species1,species2), sum, na.rm=TRUE ) 
        + tapply( shape_difference, list(species2,species1), sum, na.rm=TRUE ) )
rib.denom <- with( subset(shapediff, bone1==bone2 & bone1=="rib"), tapply( shape_difference, list(species1,species2), length )
        +  tapply( shape_difference, list(species2,species1), length  ) )
rib.speciesdiff <- rib.speciesdiff / rib.denom
pelvic.speciesdiff <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic"), tapply( shape_difference, list(species1,species2), sum, na.rm=TRUE ) 
        + tapply( shape_difference, list(species2,species1), sum, na.rm=TRUE ) )
pelvic.denom <- with( subset(shapediff, bone1==bone2 & bone1=="pelvic"), tapply( shape_difference, list(species1,species2), length )
        +  tapply( shape_difference, list(species2,species1), length  ) )
pelvic.speciesdiff <- pelvic.speciesdiff / pelvic.denom
# treedist
species.order <- match( levels(shapediff$species1), species.tree$tip.label  )
species.treedist <- treedist( species.tree )[ species.order, species.order ]

# look at this
plot( as.vector(species.treedist), as.vector(rib.speciesdiff) )
plot( as.vector(species.treedist), as.vector(pelvic.speciesdiff) )

# and resids
rib.lm <- ( lm( rib.speciesdiff[upper.tri(rib.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)] ) )
pelvic.lm <- ( lm( pelvic.speciesdiff[upper.tri(pelvic.speciesdiff)] ~ species.treedist[upper.tri(species.treedist)] ) )


####
# get values of C on internal branches
fullmat <- make.fullmat( estpar )
# posterior mean trait values relative to the root at all nodes:
fullmean <- sweep( thedata, 2, unlist(phylomeans), "-" )
fullmean[!havedata] <- fullmat[!havedata,havedata] %*% solve( fullmat[havedata,havedata], fullmean[havedata] )
fullmean <- sweep( fullmean, 2, unlist(phylomeans), "+" )
# associate mean value of adjacent nodes to each edge
edge.values <- t( apply( tree$edge, 1, function (kk) { colMeans(fullmean[kk,]) } ) )
colnames(edge.values) <- colnames(fullmean)
edge.testes <- edge.values[,"actual_testes_mass_max"]
# and testes-weighted relative time in the tree:
testes.treedist <- treedist( tree, edge.length=scale(edge.testes) )
rownames(testes.treedist) <- colnames(testes.treedist) <- rownames(thedata)
testes.spdist <- testes.treedist[match(levels(shapediff$species1),rownames(thedata)),match(levels(shapediff$species1),rownames(thedata))]

# visualization
adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
plot(adjtree,edge.width=edge.testes)
plot(adjtree,edge.width=fullmean[edge.indices])

##########
# compare residuals to testes-weighted distance
resid.rib.lm <- lm( resid(rib.lm) ~  testes.spdist[upper.tri(testes.spdist)][-rib.lm$na.action] )
resid.pelvic.lm <- lm( resid(pelvic.lm) ~  testes.spdist[upper.tri(testes.spdist)][-pelvic.lm$na.action] )
pdf( file="ribs-correlate.pdf", width=8, height=5, pointsize=10 )
layout(t(1:2))
plot( testes.spdist[upper.tri(testes.spdist)][-rib.lm$na.action], resid(rib.lm), xlab="relative testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="ribs" )
abline(coef(resid.rib.lm))
plot( testes.spdist[upper.tri(testes.spdist)][-pelvic.lm$na.action], resid(pelvic.lm), xlab="relative testes-weighted tree distance", ylab="residuals of bone distance accounting for tree distance", main="pelvic" )
abline(coef(resid.pelvic.lm))
dev.off()
### LOOKS GOOD!!
########


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
    # par = ks, sigma2S, gammaP, xi2P
    ks <- par[1]
    sigma2S <- par[2]
    gammaP <- par[3]
    xi2P <- par[4]
    return( internal.lengths * ( sigma2S + gammaP * edge.testes ) + tip.lengths * xi2P )
}

treemat <- treedist( adjtree, descendants=descendants )

# quick-and-dirty parameters
# within-species variance
tip.ribvar <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), tapply( shape_difference, species1, mean, na.rm=TRUE ) )
tip.pelvicvar <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), tapply( shape_difference, species1, mean, na.rm=TRUE ) )
# estimate of xi2P,R
xi2P <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), var( shape_difference, na.rm=TRUE ) / mean( shape_difference, na.rm=TRUE ) ) / 2
xi2R <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), var( shape_difference, na.rm=TRUE ) / mean( shape_difference, na.rm=TRUE ) ) / 2
# estimate of ksP,R
ksP <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), 2 * mean( shape_difference, na.rm=TRUE )^2 / var( shape_difference, na.rm=TRUE ) )
ksR <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), 2 * mean( shape_difference, na.rm=TRUE )^2 / var( shape_difference, na.rm=TRUE ) )
# check:
tmp <- with( subset( shapediff, species1==species2 & bone1=="pelvic" & bone2=="pelvic"), hist( shape_difference, main="pelvic", breaks=40 ) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/xi2P, df=ksP )) )
tmp <- with( subset( shapediff, species1==species2 & bone1=="rib" & bone2=="rib"), hist( shape_difference, main='rib', breaks=40 ) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/xi2R, df=ksR )) )
lines( tmp$mids, sum(tmp$counts) * diff(pchisq( tmp$breaks/xi2R*(ksP/ksR), df=ksP )), col='red' )  # looks better
# estimates
ks <- ksP
xi2P <- xi2P
xi2R <- xi2R * ksP/ksR

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

