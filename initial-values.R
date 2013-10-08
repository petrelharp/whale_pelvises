if (!exists("scriptdir")) {scriptdir <- "."}
source(paste(scriptdir,"correlated-traits-fns.R",sep="/"))
require(Matrix)

tree_file <- paste(scriptdir,"McGowenetal2009_Cetacea.modified.tree",sep='/')
species_tree<-read.nexus(file=tree_file)

bones <- read.table("62_add_centroids.out", header=TRUE)
bones <- subset(bones, ! species %in% c("ORCINUS_ORCA") )

species <- read.table("52_sexual_dimorphism.out", header=TRUE)
# cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt
morphology <- read.table("morphology_table_2013_June_27-plr.txt", sep='\t', header=TRUE)
allspecies <- sort( unique( levels(bones$species) ) )
bones$species <- factor( bones$species, levels=allspecies )
tmp <- data.frame( 
        species=allspecies,
        bodylength=tapply( bones$bodylength, bones$species, mean )
        )
stopifnot( any( tapply( bones$bodylength, bones$species, var )>0, na.rm=TRUE ) )
species <- subset( species, species %in% allspecies )
species$species <- factor( species$species, levels=allspecies )
species <- merge( species, tmp, by="species", all.x=TRUE, all.y=TRUE )
morphology$species <- factor( toupper(morphology$species), levels=allspecies )
species <- merge( species, morphology, by='species', all.x=TRUE )

# rearrange to have one entry per whale.
whichbone <- factor( paste(bones$side, bones$bone, sep='.') )
usevars <- setdiff(names(bones),c("side","bone"))
by.bone <- tapply( 1:nrow(bones), whichbone, function (k) bones[k,,drop=FALSE] )
for (k in seq_along(by.bone)) {
    by.bone[[k]] <- by.bone[[k]][usevars]
    names(by.bone[[k]])[match("centroid",names(by.bone[[k]]))] <- names(by.bone)[k]
}
whales <- by.bone[[1]]
for (k in 2:length(by.bone)) { whales <- merge( whales, by.bone[[k]], all=TRUE ) }

species_tree<-read.nexus(file=tree_file)
treedists <- dist.nodes( species_tree )
treedists <- treedists[ 1:Ntip(species_tree), 1:Ntip(species_tree) ]
rownames(treedists) <- colnames(treedists) <- species_tree$tip.label
distcats <- cut( as.vector(treedists), breaks=seq(0,.5,length.out=100) )
distpoints <- seq(0,.5,length.out=100)[-1] - diff(seq(0,.5,length.out=100))/2

####
# estimate delta

layout(1:3)
species.lengths <- tapply(bones$bodylength,bones$species,mean,na.rm=TRUE)
omit.species <- c("PHOCOENA_PHOCOENA",names(species.lengths[species.lengths>400]))
deltas <- lapply( levels(bones$bone), function (thisbone) {
    with( subset(bones,bone==thisbone & ! species %in% omit.species ), {
                plot( centroid ~ bodylength, col=species, log='xy', main=thisbone )
                thislm <- lm( log(centroid) ~ log(bodylength) + 0  )
                # abline( coef(thislm)[1]/log(10), coef(thislm)[2] )
                abline( 0, coef(thislm)[1] )
                thislm } )
} )
t.delta <- with( species, {
            plot( actual_testes_mass_max ~ bodylength, col=species, log='xy', main='testes' )
            thislm <- lm( log(actual_testes_mass_max) ~ log(bodylength) + 0 )
            abline( 0, coef(thislm)[1] )
            thislm } )

deltas
t.delta
# delta = 1.3 seems pretty good for all three...
# so sqrt(delta) = 1.14

####
# estimate sigmaL

layout(1:2)
species.lengths <- tapply(bones$bodylength,bones$species,mean,na.rm=TRUE)
length.diffs <- treedists
length.diffs[] <- NA
tmp.lengths <-  ( outer( log(species.lengths), log(species.lengths), "-" ) )
length.diffs[ rownames(tmp.lengths), colnames(tmp.lengths) ] <- tmp.lengths
length.vars <- tapply( length.diffs, distcats, var, na.rm=TRUE )
lvar.lm <- lm( (length.vars) ~ distpoints + 0 )

plot( distpoints, sqrt(length.vars), col='red', pch=20, ylim=range(c(sqrt(length.vars),length.diffs),na.rm=TRUE) )
points(treedists, length.diffs)
lines( distpoints, sqrt( coef(lvar.lm)*distpoints ), col='red' )

lvar.lm
# sigmalL = sqrt(10) = 3.16 


###
# estimate sigmaR, phylogenetic rib variance after accounting for length

rib.length.lm <- with( subset(bones,bone=='rib'), lm( log(centroid) ~ log(bodylength) ) )
rib.length.resids <- with( subset(bones,bone=='rib'), { x <- resid(rib.length.lm); names(x) <- species[-rib.length.lm$na.action]; x } )
species.rib.length.resids <- tapply( rib.length.resids, names(rib.length.resids), mean, na.rm=TRUE )
rib.resid.diffs <- treedists
rib.resid.diffs[] <- NA
tmp <- outer( species.rib.length.resids, species.rib.length.resids, "-" )
rib.resid.diffs[ rownames(tmp), colnames(tmp) ] <- tmp 
rib.resid.vars <- tapply( rib.resid.diffs, distcats, var, na.rm=TRUE )
rib.resid.lm <- lm( (rib.resid.vars) ~ distpoints + 0 )

plot( treedists, rib.resid.diffs  )
points( distpoints, sqrt( rib.resid.vars ), col='red', pch=20 )
lines( distpoints, sqrt(coef(rib.resid.lm)*distpoints), col='red' )

rib.resid.lm

# less phylogenetic signal
# sigmaR = sqrt(.2478)) = .5

###
# estimate sigmaP, phylogenetic pelvic variance after accounting for length

pelvic.length.lm <- with( subset(bones,bone=='pelvic'), lm( log(centroid) ~ log(bodylength) ) )
pelvic.length.resids <- with( subset(bones,bone=='pelvic'), { x <- resid(pelvic.length.lm); names(x) <- species[-pelvic.length.lm$na.action]; x } )
species.pelvic.length.resids <- tapply( pelvic.length.resids, names(pelvic.length.resids), mean, na.rm=TRUE )
pelvic.resid.diffs <- treedists
pelvic.resid.diffs[] <- NA
tmp <- outer( species.pelvic.length.resids, species.pelvic.length.resids, "-" )
pelvic.resid.diffs[ rownames(tmp), colnames(tmp) ] <- tmp 
pelvic.resid.vars <- tapply( pelvic.resid.diffs, distcats, var, na.rm=TRUE )
pelvic.resid.lm <- lm( (pelvic.resid.vars) ~ distpoints + 0 )

plot( treedists, pelvic.resid.diffs  )
points( distpoints, sqrt( pelvic.resid.vars ), col='red', pch=20 )
lines( distpoints, sqrt(coef(pelvic.resid.lm)*distpoints), col='red' )

pelvic.resid.lm

# more phylogenetic signal
# sigmaP = sqrt(1.26) = 1.12

###
# and betaT, phylogenetic testes size after accounting for length

testes.length.lm <- with( species, lm( log(actual_testes_mass_max) ~ log(bodylength) ) )
testes.length.resids <- with( species, { x <- resid(testes.length.lm); names(x) <- species[-testes.length.lm$na.action]; x } )
testes.resid.diffs <- treedists
testes.resid.diffs[] <- NA
tmp <- outer( testes.length.resids, testes.length.resids, "-" )
testes.resid.diffs[ rownames(tmp), colnames(tmp) ] <- tmp 
testes.resid.vars <- tapply( testes.resid.diffs, distcats, var, na.rm=TRUE )
testes.resid.lm <- lm( (testes.resid.vars) ~ distpoints + 0 )

plot( treedists, testes.resid.diffs  )
points( distpoints, sqrt( testes.resid.vars ), col='red', pch=20 )
lines( distpoints, sqrt(coef(testes.resid.lm)*distpoints), col='red' )

testes.resid.lm

# YES phylogenetic signal
# betaT = sqrt(41.65) = 6.5


###
# and betaP, something like covariance between testes resids and pelvic resids?

plot( testes.length.resids, species.pelvic.length.resids[names(testes.length.resids)] )
identify( testes.length.resids, species.pelvic.length.resids[names(testes.length.resids)], labels=names(testes.length.resids) )
lm( species.pelvic.length.resids[names(testes.length.resids)] ~ testes.length.resids  )
# correlated, coefficient 0.12
#  driven most by c("EUBALAENA_GLACIALIS","PONTOPORIA_BLAINVILLEI","MESOPLODON_CARLHUBBSI") ??
#   and also by c("ESCHRICHTIUS_ROBUSTUS","PHOCOENA_PHOCOENA", "ZIPHIUS_CAVIROSTRIS", "BALAENOPTERA_PHYSALUS")


###
# what about testes-rib residual correlation?
plot( testes.length.resids, species.rib.length.resids[names(testes.length.resids)] )
identify( testes.length.resids, species.rib.length.resids[names(testes.length.resids)], labels=names(testes.length.resids) )
lm( species.rib.length.resids[names(testes.length.resids)] ~ testes.length.resids  )
# NOT correlated.


####
# estimate zetaL: within species, length
length.resids <- log(bones$bodylength) - log(species.lengths[ match(bones$species,names(species.lengths)) ])
sqrt(var(length.resids))
plot(bones$bodylength, length.resids,log='x')

# zetaL = .05

####
# estimate zetaR: within species, rib
species.ribs <- with(subset(bones,bone=='rib'), tapply( centroid, species, mean, na.rm=TRUE ) )
rib.resids <- with(subset(bones,bone=='rib'), log(centroid) - log(species.ribs[ match(species,names(species.ribs)) ]) )
sqrt(var(rib.resids,na.rm=TRUE))
with(subset(bones,bone=='rib'), plot(centroid, rib.resids,log='x') )

# zetaR = .06

###
# and omegaR, within indivs, rib
indiv.ribs <- with(subset(bones,bone=='rib'), tapply( centroid, specimen, mean, na.rm=TRUE ) )
rib.indiv.resids <- with(subset(bones,bone=='rib'), log(centroid) - log(indiv.ribs[ match(specimen,names(indiv.ribs)) ]) )
sqrt(var(rib.indiv.resids,na.rm=TRUE))
with(subset(bones,bone=='rib'), plot( centroid, rib.indiv.resids, log='x' ) )

# omegaR = .012

####
# and zetaP: pelvis
species.pelvics <- with(subset(bones,bone=='pelvic'), tapply( centroid, species, mean, na.rm=TRUE ) )
pelvic.resids <- with(subset(bones,bone=='pelvic'), log(centroid) - log(species.pelvics[ match(species,names(species.pelvics)) ]) )
sqrt(var(pelvic.resids,na.rm=TRUE))
with(subset(bones,bone=='pelvic'), plot(centroid, pelvic.resids,log='x') )

# zetaP = .12

###
# and omegaP, within indivs, pelvic
indiv.pelvics <- with(subset(bones,bone=='pelvic'), tapply( centroid, specimen, mean, na.rm=TRUE ) )
pelvic.indiv.resids <- with(subset(bones,bone=='pelvic'), log(centroid) - log(indiv.pelvics[ match(specimen,names(indiv.pelvics)) ]) )
sqrt(var(pelvic.indiv.resids,na.rm=TRUE))
with(subset(bones,bone=='pelvic'), plot( centroid, pelvic.indiv.resids, log='x' ) )

# omegaP = .02
