#/usr/bin/R --vanilla


## simulation parameters:
# dimensionality
ks <- nvars <- 20
# within species and individual variance
# xi2s <- 0.1
# xi2i <- 0.05
xi2s <- 0.001
xi2i <- 0.0005

require(ape)
require(Matrix)

basedir <- scriptdir <- ".."
source(paste(scriptdir,"correlated-traits-fns.R",sep="/"))

tree_file <- paste(basedir, "consensusTree_ALL_CETACEA.tree", sep='/')
species_tree<-read.nexus(file=tree_file)
load( paste(basedir,"all-sample-tree.RData",sep='/') ) # gets tree

whales <- read.csv(paste(basedir,"whales.csv",sep='/'),header=TRUE)
allspecies <- intersect(tree$tip.label, species_tree$tip.label)

do.simple <- FALSE
if (do.simple) {
    ## SIMPLE EXAMPLE
    allspecies <- c( "PHOCOENOIDES_DALLI", "STENELLA_ATTENUATA", "STENELLA_LONGIROSTRIS" )
    whales <- droplevels( subset(whales,species%in%allspecies) )
    tree <- drop.tip( tree, setdiff( tree$tip.label, c( allspecies, as.character(whales$specimen) ) ) )

    bones <- read.table(paste(basedir,"62_add_centroids.out",sep='/'), header=TRUE)
    bones <- droplevels( subset(bones, species %in% allspecies & ! (specimen == "LACM_54109") ) )
    bones <- bones[ setdiff( colnames(bones), "absolute_volume" ) ]
    species <- read.table(paste(basedir,"52_sexual_dimorphism.out",sep='/'), header=TRUE)
    morphology <- read.table(paste(basedir,"morphology_table_2013_June_27-plr.txt",sep='/'), sep='\t', header=TRUE)
    stopifnot( all.equal( allspecies , sort( unique( levels(bones$species) ) ) ) )
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

    # get thedata ready
    thedata <- matrix( c(NA), nrow=Nnode(tree)+Ntip(tree), ncol=6 )
    colnames(thedata) <- c("bodylength","actual_testes_mass_max","left.rib","right.rib","left.pelvic","right.pelvic")
    rownames(thedata) <- c( tree$tip.label, paste("NA",1:Nnode(tree),sep='.') )
    species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
    sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
    thedata[species.nodes,intersect(names(species),colnames(thedata))] <- as.matrix( species[match(rownames(thedata)[species.nodes],species$species),intersect(names(species),colnames(thedata))] )
    thedata[sample.nodes,intersect(names(whales),colnames(thedata))] <- as.matrix( whales[match(rownames(thedata)[sample.nodes],whales$specimen),intersect(names(whales),colnames(thedata))] )
    thedata <- log(thedata)
    havedata <- !is.na(thedata)
} else {
    load(paste(basedir,"thedata-and-covmatrices.Rdata",sep='/'))
}

##
# the actual data
shapediff <- read.table("../61_pairwise_shapes.out", header=TRUE)
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



# from parse-correlated-data.R
## tree stuff:
descendants <- get.descendants(tree)  # [i,j] TRUE means that j is below i
n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
rootnode <- which.max(n.offspring)
# edge.indices[k] is the node associated with edge k
# match(j,edge.indices) is the edge associated with node j
edge.indices <- tree$edge[,2] # associate each edge with the downstream node
tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
internal.edges <- setdiff( 1:Nedge(tree), tip.edges )
species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
internal.nodes <- Ntip(tree) + (1:Nnode(tree))


# [k,j] TRUE if node j is below edge k
edge.descendants <- descendants[edge.indices,]

# normalizing whatnot
adjtree <- tree
adjtree$edge.length[tip.edges] <- xi2s

fake.noise <- function () {
    edgediffs <- rnorm( nvars * Nedge(adjtree) )
    dim(edgediffs) <- c(nvars,Nedge(adjtree))
    edgediffs <- sweep( edgediffs, 2, sqrt(adjtree$edge.length), "*" )
    return(edgediffs)
}

###
fake.shapes <- function () {
    edgediffs <- fake.noise()
    # sum down the tree, getting shapes for each sample
    # rows correspond to nodes in the tree
    simdata <- t( edgediffs %*% edge.descendants )
    lr.noise <- sqrt(xi2i) * rnorm( 2 * nrow(simdata) * nvars )
    dim(lr.noise) <- c( 2 * nrow(simdata), nvars )
    simdata <- simdata[ rep(1:nrow(simdata),each=2), ] + lr.noise
    rownames(simdata) <- paste( rep(c( adjtree$tip.label, paste("node",1:Nnode(tree),sep='') ), each=2), rep(c("L","R"),nrow(simdata)/2), sep="." )
    return(simdata)
}

###
fake.data <- function () {
    simdata <- fake.shapes()
    diffs <- rowSums( apply( simdata, 2, function (x) outer(x,x,"-") )^2 )
    dim(diffs) <- c(nrow(simdata),nrow(simdata))
    spsides <- strsplit( rownames(simdata), ".", fixed=TRUE )
    specimens <- sapply( spsides, "[[", 1 )
    sides <- sapply( spsides, "[[", 2 )
    ut <- upper.tri(diffs)
    diffdata <- data.frame(
                specimen1 = specimens[row(diffs)[ut]],
                specimen2 = specimens[col(diffs)[ut]],
                side1 = sides[row(diffs)[ut]],
                side2 = sides[col(diffs)[ut]],
                sqdist = diffs[ut]
            )
    return(diffdata)
}

simdata <- fake.data()

write.csv(simdata,file='full-shape-simdata.csv',row.names=FALSE)

## many datasets (2m for 1000)
simdatas <- replicate(1000, { simdata <- fake.shapes(); diffs <- rowSums( apply( simdata, 2, function (x) outer(x,x,"-") )^2 ); dim(diffs) <- c(nrow(simdata),nrow(simdata)); diffs[upper.tri(diffs)]  })
save(simdata,simdatas,file="many-full-shape-simdata.RData")
