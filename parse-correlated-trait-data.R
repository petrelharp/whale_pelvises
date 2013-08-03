#!/usr/bin/R

require(ape)
source("correlated-traits-fns.R")

tree_file <- "consensusTree_ALL_CETACEA.tree"
species_tree<-read.nexus(file=tree_file)

bones <- read.table("62_add_centroids.out", header=TRUE)
bones <- droplevels( subset(bones, ! species %in% c("ORCINUS_ORCA") ) )
bones <- bones[ setdiff( colnames(bones), "absolute_volume" ) ]

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

if (FALSE) { ## for testing
    whales <- droplevels( subset( whales, species %in% c( "BALAENOPTERA_PHYSALUS", "BALAENOPTERA_MUSCULUS" ) ) )
}

# set up the tree
within_length <- 1
tree <- species_tree
for (sp in intersect(levels(whales$species),tree$tip.label)) {
    theseones <- (whales$species==sp)
    nsamps <- sum(theseones)
    tmptree <- list(  # include zero-length edge (to first node) for species observation
        edge = cbind( rep(nsamps+2,nsamps+1), 1:(nsamps+1) ),
        edge.length = c(0,rep(within_length,nsamps)),
        Nnode = 1,
        tip.label = c(sp,levels(whales$specimen)[as.numeric(whales$specimen[theseones])]),
        root.edge = 0  # this is not a length.
        )
    class(tmptree) <- "phylo"
    tree <- bind.tree(x=tree,y=tmptree,where=match(sp,tree$tip.label))
}

tree <- drop.tip( tree, setdiff( tree$tip.label, c(levels(whales$specimen),levels(whales$species)) ) )
stopifnot( setequal( c(levels(whales$specimen),levels(whales$species)), tree$tip.label ) )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% 1:Ntip(tree) ] %in% c(0,within_length) ) )

## tree stuff:
descendants <- get.descendants(tree)
n.offspring <- rowSums(descendants)  # number of tips each edge contributes to
rootnode <- which.max(n.offspring)
# edge.indices[k] is the node associated with edge k
# match(j,edge.indices) is the edge associated with node j
edge.indices <- tree$edge[,2] # associate each edge with the downstream node
tip.edges <- match( 1:Ntip(tree), edge.indices )  # which edges correspond to tips... note these are the first Ntip(tree) nodes
species.nodes <- match( levels(whales$species), tree$tip.label )  # which nodes are zero-length "species" nodes
sample.nodes <- setdiff( 1:Ntip(tree), species.nodes )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% species.nodes ] == 0 ) )
stopifnot( all( tree$edge.length[ tree$edge[,2] %in% sample.nodes ] == within_length ) )
stopifnot( ! ( rootnode %in% edge.indices ) )

###
# Get the data all together: [i,j] is j-th variable for i-th node in the tree
# again: the variables are, in order: Length, Testes, Rib-left, Rib-right, Pelvis-left, Pelvis-right
thedata <- matrix( c(NA), nrow=Nnode(tree)+Ntip(tree), ncol=6 )
colnames(thedata) <- c("bodylength","actual_testes_mass_max","left.rib","right.rib","left.pelvic","right.pelvic")
rownames(thedata) <- c( tree$tip.label, paste("NA",1:Nnode(tree),sep='.') )
# species obs
thedata[species.nodes,intersect(names(species),colnames(thedata))] <- as.matrix( species[match(rownames(thedata)[species.nodes],species$species),intersect(names(species),colnames(thedata))] )
thedata[sample.nodes,intersect(names(whales),colnames(thedata))] <- as.matrix( whales[match(rownames(thedata)[sample.nodes],whales$specimen),intersect(names(whales),colnames(thedata))] )
thedata <- log(thedata)
havedata <- !is.na(thedata)

# "phylogenetic" mean: use imputed data to deal with missing values
imputed.data <- thedata
for (sp in levels(whales$species)) {
    is.sp <- ( rownames(imputed.data) == sp ) | ( whales$species[ match( rownames(imputed.data), whales$specimen ) ] == sp )
    is.sp[ is.na(is.sp) ] <- FALSE
    sp.means <- colMeans( imputed.data[is.sp,], na.rm=TRUE )
    if (any(is.na(sp.means))) {
        # NA-out those species with some variables still missing
        imputed.data[is.sp] <- NA
    } else {
        imputed.data[is.sp & is.na(imputed.data)] <- sp.means[ col( imputed.data )[is.sp & is.na(imputed.data)] ]
    }
}


adjtree <- tree
adjtree$edge.length[tip.edges] <- (.05/3.16)^2  # reasonable value from initial-values.R
phylomeans <- lapply( 1:ncol(imputed.data), function(k) phylomean(imputed.data[1:Ntip(tree),k], tree=adjtree) )
names(phylomeans) <- colnames(imputed.data)
edgeweights <- tip.to.edge.weights( attr(phylomeans[[1]],"weights"), adjtree, descendants )

# and, the weights of each species in the phylogenetic mean:
# normalize by "phylogenetic" mean
thedata <- sweep( thedata, 2, unlist(phylomeans), "-" )

if (FALSE) {  # DO THIS LATER
    # normalize by sexual dimorphism
    data.specimens <- match( rownames(thedata), whales$specimen )
    females <- ( whales$sex[ data.specimens ] == "F" )
    data.species <- whales$species[match(data.specimens,whales$specimen)]
    renorms <- species$sexual_size_dimorphism[match(data.species,species$species)]
    for (x in c("bodylength",levels(whichbone))) {
        whales[[x]][ whales$sex=="F" ] <- whales[[x]][ whales$sex=="F" ] / renorms[ whales$sex=="F" ]
    }
}


#####
# covariance matrices:

# "internal" edges --
internal.lengths <- tree$edge.length
internal.lengths[ tip.edges ] <- 0
# lengths of shared internal branches leading to the root for all pairs of nodes in the tree
species.shared.brlens <- shared.branchlengths( tree, internal.lengths*(1-2*edgeweights), descendants )
# distances in internal branches in the tree for all pairs of nodes
species.treedist <- treedist( tree, internal.lengths*edgeweights, descendants )
species.treemat <- ( species.shared.brlens - species.treedist + sum(edgeweights^2 * internal.lengths) )

# and, in the tips
tip.lengths <- tree$edge.length
tip.lengths[ - tip.edges ] <- 0
tip.treedist <- treedist( tree, tip.lengths*edgeweights, descendants )
sample.treemat <- ( 0 - tip.treedist + sum(edgeweights^2 * tip.lengths) )


###
# write out
save( species.treemat, sample.treemat, thedata, file="thedata-and-covmatrices.Rdata" )

# write.csv( treedist, file="all-sample-treedist.csv", row.names=FALSE)
# write.csv( tipdist, file="all-sample-tipdist.csv", row.names=FALSE)
# write.csv( thedata, file="all-data-rejiggered.csv", row.names=TRUE)
# write.csv( phylomeans, file="phylomeans.csv", row.names=FALSE )
# write.tree( tree, file="all-sample-tree.R")
