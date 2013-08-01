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

#####
## ok do a big linear model
require(plyr)

whmeans <- ddply( whales, "species", summarise, bodylength=mean(bodylength,na.rm=TRUE), left.pelvic=mean(left.pelvic,na.rm=TRUE), left.rib=mean(left.rib,na.rm=TRUE), right.pelvic=mean(right.pelvic,na.rm=TRUE), right.rib=mean(right.rib,na.rm=TRUE) )
whsds <- ddply( whales, "species", summarise, bodylength=sqrt(var(bodylength,na.rm=TRUE)), left.pelvic=sqrt(var(left.pelvic,na.rm=TRUE)), left.rib=sqrt(var(left.rib,na.rm=TRUE)), right.pelvic=sqrt(var(right.pelvic,na.rm=TRUE)), right.rib=sqrt(var(right.rib,na.rm=TRUE)) )

layout(matrix(1:4,nrow=2))
with(whmeans, plot( left.pelvic, right.pelvic ) ); abline(0,1)
with(whmeans, plot( left.rib, right.rib ) ); abline(0,1)

##
lm0 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + sex ) )
lm1 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + species ) )
lm2 <- with( bones, lm( log(absolute_volume) ~ log(bodylength) + species + sex ) )


