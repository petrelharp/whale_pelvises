source("correlated-traits-fns.R")
require(Matrix)

tree_file <- "consensusTree_ALL_CETACEA.tree"
species_tree<-read.nexus(file=tree_file)

bones <- read.table("50_make_datamatrix.out", header=TRUE)
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
    names(by.bone[[k]])[match("absolute_volume",names(by.bone[[k]]))] <- names(by.bone)[k]
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
                plot( (absolute_volume) ~ (bodylength), col=species, log='xy', main=thisbone )
                thislm <- lm( log(absolute_volume) ~ log(bodylength) + 0  )
                # abline( coef(thislm)[1]/log(10), coef(thislm)[2] )
                abline( 0, coef(thislm)[1] )
                thislm } )
} )
t.delta <- with( species, {
            plot( actual_testes_mass_max ~ bodylength, col=species, log='xy', main='testes' )
            thislm <- lm( log(actual_testes_mass_max) ~ log(bodylength) + 0 )
            abline( 0, coef(thislm)[1] )
            thislm } )
# delta = 1.75 seems pretty good for all three...
# so sqrt(delta) = 1.32

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

rib.length.lm <- with( subset(bones,bone=='rib'), lm( log(absolute_volume) ~ log(bodylength) ) )
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

# not much phylogenetic signal
# sigmaR = sqrt(.414) = .64

###
# estimate sigmaP, phylogenetic pelvic variance after accounting for length

pelvic.length.lm <- with( subset(bones,bone=='pelvic'), lm( log(absolute_volume) ~ log(bodylength) ) )
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

# YES phylogenetic signal
# sigmaP = sqrt(5.3) = 2.3



####
# estimate zetaL: within species, length
length.resids <- log(bones$bodylength) - log(species.lengths[ match(bones$species,names(species.lengths)) ])
sqrt(var(length.resids))
plot(bones$bodylength, length.resids,log='x')

# zetaL = .05

####
# estimate zetaR: within species, rib
species.ribs <- with(subset(bones,bone=='rib'), tapply( absolute_volume, species, mean, na.rm=TRUE ) )
rib.resids <- with(subset(bones,bone=='rib'), log(absolute_volume) - log(species.ribs[ match(species,names(species.ribs)) ]) )
sqrt(var(rib.resids,na.rm=TRUE))
with(subset(bones,bone=='rib'), plot(absolute_volume, rib.resids,log='x') )

# zetaR = .144

###
# and omegaR, within indivs, rib
indiv.ribs <- with(subset(bones,bone=='rib'), tapply( absolute_volume, specimen, mean, na.rm=TRUE ) )
rib.indiv.resids <- with(subset(bones,bone=='rib'), log(absolute_volume) - log(indiv.ribs[ match(specimen,names(indiv.ribs)) ]) )
sqrt(var(rib.indiv.resids,na.rm=TRUE))
with(subset(bones,bone=='rib'), plot( absolute_volume, rib.indiv.resids, log='x' ) )

# omegaR = .036

####
# and zetaP: pelvis
species.pelvics <- with(subset(bones,bone=='pelvic'), tapply( absolute_volume, species, mean, na.rm=TRUE ) )
pelvic.resids <- with(subset(bones,bone=='pelvic'), log(absolute_volume) - log(species.pelvics[ match(species,names(species.pelvics)) ]) )
sqrt(var(pelvic.resids,na.rm=TRUE))
with(subset(bones,bone=='pelvic'), plot(absolute_volume, pelvic.resids,log='x') )

# zetaP = .29

###
# and omegaR, within indivs, rib
indiv.pelvics <- with(subset(bones,bone=='pelvic'), tapply( absolute_volume, specimen, mean, na.rm=TRUE ) )
pelvic.indiv.resids <- with(subset(bones,bone=='pelvic'), log(absolute_volume) - log(indiv.pelvics[ match(specimen,names(indiv.pelvics)) ]) )
sqrt(var(pelvic.indiv.resids,na.rm=TRUE))
with(subset(bones,bone=='pelvic'), plot( absolute_volume, pelvic.indiv.resids, log='x' ) )

# omegaR = .035

######


lms <- with( subset(bones,bone=="pelvic"), list( 
                lm( log(absolute_volume) ~ log(bodylength) ),
                lm( log(absolute_volume) ~ specimen ) )
        )

#####
## ok do a big linear model
require(plyr)

whmeans <- ddply( whales, "species", summarise, bodylength=mean(bodylength,na.rm=TRUE), left.pelvic=mean(left.pelvic,na.rm=TRUE), left.rib=mean(left.rib,na.rm=TRUE), right.pelvic=mean(right.pelvic,na.rm=TRUE), right.rib=mean(right.rib,na.rm=TRUE) )
whsds <- ddply( whales, "species", summarise, bodylength=sqrt(var(bodylength,na.rm=TRUE)), left.pelvic=sqrt(var(left.pelvic,na.rm=TRUE)), left.rib=sqrt(var(left.rib,na.rm=TRUE)), right.pelvic=sqrt(var(right.pelvic,na.rm=TRUE)), right.rib=sqrt(var(right.rib,na.rm=TRUE)) )

layout(matrix(1:4,nrow=2))
with(whmeans, plot( left.pelvic, right.pelvic ) ); abline(0,1)
with(whmeans, plot( left.rib, right.rib ) ); abline(0,1)

lm0 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + sex ) )
lm1 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + species ) )
lm2 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + species + sex ) )


