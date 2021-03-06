#!/usr/bin/R

bones <- read.table("62_add_centroids.out", header=TRUE)
shapediff <- read.table("61_pairwise_shapes.out", header=TRUE)
species <- read.table("52_sexual_dimorphism.out", header=TRUE)
# cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt
morphology <- read.table("morphology_table_2013_June_27-plr.txt", sep='\t', header=TRUE)

# EXCLUDE PROBLEM BONE:
problems <- list( list(specimen='USNM_504107', bone='rib', side='right') )
for (x in problems) {
    bones <- subset( bones, ! ( (specimen == x$specimen) & (bone == x$bone) & (side == x$side) ) )
}
tmp.bones <- with( bones, paste( specimen, bone, side ) )
tmp.shape1 <- with( shapediff, paste( specimen1, bone1, side1 ) )
tmp.shape2 <- with( shapediff, paste( specimen2, bone2, side2 ) )
# setdiff( tmp.bones, c(tmp.shape1,tmp.shape2) )
grep( "REP", setdiff( c(tmp.shape1,tmp.shape2), tmp.bones ), invert=TRUE, value=TRUE )
shapediff <- subset( shapediff, ( tmp.shape1 %in% tmp.bones) & ( tmp.shape2 %in% tmp.bones) )


allspecies <- sort( unique( c( levels(bones$species), levels(species$species) ) ) )
bones$species <- factor( bones$species, levels=allspecies )
tmp <- data.frame( 
        species=allspecies,
        bodylength=tapply( bones$bodylength, bones$species, mean )
        )
stopifnot( any( tapply( bones$bodylength, bones$species, var )>0, na.rm=TRUE ) )
species$species <- factor( species$species, levels=allspecies )
species <- merge( species, tmp, by="species", all.x=TRUE, all.y=TRUE )

morphology$species <- factor( toupper(morphology$species), levels=allspecies )
species <- merge( species, morphology, by='species', all.x=TRUE )

specimens <- sort(levels(bones$specimen))
bones$specimen <- factor( bones$specimen, levels=specimens )
shapediff$specimen1 <- factor( shapediff$specimen1, levels=specimens )
shapediff$specimen2 <- factor( shapediff$specimen2, levels=specimens )

require(ape)
# tree <- read.nexus("consensusTree_10_species.txt")

bones$genus <- bones$species
levels(bones$genus) <- sapply( strsplit( levels(bones$genus), "_" ), "[[", 1 )

cetacean_tree <- read.nexus(file="McGowenetal2009_Cetacea.modified.tree")
tree<-drop.tip(cetacean_tree, setdiff( cetacean_tree$tip.label, allspecies ) )
species_dist <- cophenetic(tree)
species_dist <- species_dist[ levels(bones$species),levels(bones$species)]

########
## look at shape differences as a function of phylogenetic distance


shapediff$species1 <- bones$species[ match(shapediff$specimen1,bones$specimen) ] 
shapediff$species2 <- bones$species[ match(shapediff$specimen2,bones$specimen) ] 
shapediff$genus1 <- bones$genus[ match(shapediff$specimen1,bones$specimen) ]
shapediff$genus2 <- bones$genus[ match(shapediff$specimen2,bones$specimen) ]
shapediff$sex1 <- bones$sex[ match(shapediff$specimen1,bones$specimen) ]
shapediff$sex2 <- bones$sex[ match(shapediff$specimen2,bones$specimen) ]

shapediff$comparison <- factor( "wider", levels=c("intraindiv","intrasex","intraspecies","intragenus","wider") )
shapediff$comparison[ ( shapediff$genus1 == shapediff$genus2  ) ] <- "intragenus"
shapediff$comparison[ ( shapediff$species1 == shapediff$species2  ) ] <- "intraspecies"
shapediff$comparison[ (shapediff$comparison=="intraspecies") & ( shapediff$sex1 == shapediff$sex2  ) ] <- "intrasex"
shapediff$comparison[ ( shapediff$specimen1 == shapediff$specimen2 ) ] <- "intraindiv"

shapediff$phylodist <- species_dist[ cbind( shapediff$species1, shapediff$species2 ) ]

genuscolors <- adjustcolor( rainbow(nlevels(bones$genus)+5)[1:nlevels(bones$genus)], .5 )
genusadj <- (-.5) + (1/(nlevels(bones$genus)+4))*as.numeric(shapediff$genus1)
genusadj2 <- (-.5) + (1/(nlevels(bones$genus)+4))*as.numeric(shapediff$genus2)

layout(t(1:2))
for (thisbone in c("pelvic","rib")) {
    with( subset( shapediff, (bone1==bone2) & (bone1==thisbone) & (species1==species2) ), {
                # plot( as.numeric(species1)+(as.numeric(comparison)-2.5)/10, shape_difference, col=comparison )
                dotchart( shape_difference, groups=paste(as.character(species1),as.character(comparison)) )
            } )
}


pdf(file="figures/species-genus-comparisons.pdf",width=10,height=8,pointsize=10)
layout(matrix(c(1,2,3,3),nrow=2),widths=c(5,1))
for (thisbone in c("pelvic","rib")) {
    plot( shape_difference ~ comparison, data=shapediff, subset=(bone1==bone2 & bone1==thisbone), main=thisbone )
    # points( shape_difference ~ jitter(as.numeric(comparison),factor=.4), data=shapediff, subset=(comparison=="intragenus"), col=genus1 )
    # points( shape_difference ~ jitter(as.numeric(comparison),factor=.4), data=shapediff, subset=(comparison%in%c("intraspecies","intraindiv")), col=species1 )
    points( shape_difference ~ I(genusadj + as.numeric(comparison)), pch=20, col=genuscolors[genus1], data=shapediff, subset=(bone1==bone2 & bone1==thisbone) )
    points( shape_difference ~ I(genusadj2 + as.numeric(comparison)), pch=20, col=genuscolors[genus2], data=shapediff, subset=(bone1==bone2 & bone1==thisbone & comparison=="wider") )
    plot( shape_difference ~ comparison, data=shapediff, subset=(bone1==bone2 & bone1==thisbone), main=thisbone, add=TRUE )
}
opar <- par(mar=c(5,0,4,0)+.1)
plot(0,type='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend("topleft",legend=levels(bones$genus),col=genuscolors,pch=20)
par(opar)
dev.off()

pdf(file="figures/species-genus-side-comparisons.pdf",width=10,height=8,pointsize=10)
layout(matrix(c(1,2,3,3),nrow=2),widths=c(5,1))
for (thisbone in c("pelvic","rib")) {
    plot( shape_difference ~ comparison, data=shapediff, subset=(bone1==bone2 & bone1==thisbone), main=thisbone )
    points( shape_difference ~ I(genusadj + as.numeric(comparison)), pch=20, col=adjustcolor(c("red","black"),.4)[1+(side1==side2)], data=shapediff, subset=(bone1==bone2 & bone1==thisbone) )
    plot( shape_difference ~ comparison, data=shapediff, subset=(bone1==bone2 & bone1==thisbone), main=thisbone, add=TRUE )
}
opar <- par(mar=c(5,0,4,0)+.1)
plot(0,type='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend("topleft",legend=levels(bones$genus),col=genuscolors,pch=20)
legend('bottomleft', legend=c('different sides','same side'), col=c('red','black'), pch=20 )
par(opar)
dev.off()

## identify problem bones
if (FALSE) {
    shapediff$genusadj <- genusadj

    layout(t(1:2))
    for (thisbone in c("rib","pelvic")) {
        problems <- c("LACM_97405","USNM_572775")

        plot( shape_difference ~ I(genusadj + as.numeric(comparison)), pch=20, col=adjustcolor(c("black",'red'),.4)[1+(specimen1%in%problems | specimen2%in%problems)], data=shapediff, subset=(bone1==bone2 & bone1==thisbone) )
        with( subset(shapediff, (bone1==bone2 & bone1==thisbone)), identify( genusadj + as.numeric(comparison), shape_difference, labels=paste(specimen1,side1,specimen2,side2) ) )
    }

    # without problem bones
    plot( shape_difference ~ I(genusadj + as.numeric(comparison)), pch=20, col=adjustcolor(c("red","black"),.4)[1+(side1==side2)], data=shapediff, subset=(bone1==bone2 & bone1==thisbone & ! ( specimen1 %in% problems | specimen2 %in% problems ) ) )

    # who tends to have the most diffs with others?
    sort( with( subset(shapediff,comparison=="intraspecific" & specimen1 != specimen2), 
        tapply( shape_difference, specimen1, sum, na.rm=TRUE ) + tapply( shape_difference, specimen2, sum, na.rm=TRUE ) ) )
    sort( with( subset(shapediff,comparison=="intragenus" & specimen1 != specimen2), 
        tapply( shape_difference, specimen1, sum, na.rm=TRUE ) + tapply( shape_difference, specimen2, sum, na.rm=TRUE ) ) )

}


pdf(file="figures/shapediff-by-phylodist.pdf",width=10,height=8,pointsize=10)

layout(matrix(c(1,2,3,3),nrow=2),widths=c(5,1))
for (thisbone in c("pelvic","rib")) {
    with( subset(shapediff,bone1==bone2&bone1==thisbone), plot( shape_difference ~ phylodist, col=genuscolors[genus1], main=thisbone, xlim=c(0,75) ) )
}
opar <- par(mar=c(5,0,4,0)+.1)
plot(0,type='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend("topleft",legend=levels(bones$genus),col=genuscolors,pch=20)
par(opar)

layout(matrix(c(1,2,3,3),nrow=2),widths=c(5,1))
for (thisbone in c("pelvic","rib")) {
    with( subset(shapediff,bone1==bone2&bone1==thisbone), plot( shape_difference ~ phylodist, col=genuscolors[genus1], main=paste(thisbone,"zoomed in"), xlim=c(0,20) ) )
}
opar <- par(mar=c(5,0,4,0)+.1)
plot(0,type='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend("topleft",legend=levels(bones$genus),col=genuscolors,pch=20)
par(opar)

dev.off()

